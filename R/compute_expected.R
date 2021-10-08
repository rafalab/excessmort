#' Compute expected counts for each day
#' 
#' Compute the expected death count for each unit of time. 
#' We assume counts are ovedsispersed Poisson distributed with a 
#' trend that accounts for changes in death rate across time and a seasonal effect. The function take data frame with 
#' dates and counts and returns the data frame with the expected counts as a new 
#' column. It also returns a logical column that is `TRUE` if that entry was 
#' used in the estimation procedure.
#' 
#' @param counts A data frame with dates, counts, and population size
#' @param exclude A list of dates to exclude when fitting the model
#' @param trend.knots.per.year Number of knots per year used for the time trend
#' @param harmonics Number of harmonics to include in the seasonal effect
#' @param frequency Number of data points per year. If not provided, the function attempts to estimate it
#' @param weekday.effect A logical that determines if a day of the week effect is included in the model
#' @param keep.components A logical that if `TRUE` forces the function to return the estimated trend, seasonal model, and weekday effect, if included in the model.
#' @param verbose A logical that if `TRUE` makes function prints out updates on the estimation procedure
#' 
#' @return The `counts` data.frame with two columns added: `expected` and `excluded`. 
#' The `expected` column is the estimated expected value of the counts for that date.
#' The `excluded` column is a logical vector denoting if that date was excluded when
#' estimating the expected value.
#' 
#' If the argument `keep.components` is `TRUE` a list is returned with `counts`
#' data.frame in the first component, the estimated trend in the second, the 
#' estimated seasonal effect in the third and the estimated weekday effects in the fourth.
#' 
#' 
#' @examples
#' data(new_jersey_counts)
#' exclude_dates <- as.Date("2012-10-29") + 0:180
#' counts <- compute_expected(new_jersey_counts, exclude = exclude_dates, weekday.effect = TRUE)
#' library(ggplot2)
#' expected_plot(counts)
#' 
#' @export
#' @importFrom stats glm contr.sum model.matrix contrasts<-
#'
compute_expected <- function(counts,
                             exclude = NULL,
                             trend.knots.per.year = 1/7,
                             harmonics = 2,
                             frequency = NULL,
                             weekday.effect = FALSE,
                             keep.components = FALSE,
                             verbose = TRUE){

  ## check column names
  if(!all(c("date", "outcome", "population") %in% names(counts))) stop("counts must have columns named date, outcome, and poulation.")

  if(any(table(counts$date))>1) stop("Each date can appear at most once.")

  if(any(c("expeceted", "excluded") %in% names(counts))) warning("expected and excluded columns will be overwritten.")

  if(any(is.na(counts$date)) | any(is.na(counts$outcome)) | any(is.na(counts$population)))
    stop("No NAs permited in date, outcome, or population columns.")

  if(!lubridate::is.Date(counts$date)) stop("date column must be class Date.")

  if(!is.numeric(counts$outcome)) stop("outome column must be counts.")

  if(!is.numeric(counts$population)) stop("population column must be numeric.")

  if(is.null(exclude)){
    warning("No dates excluded. We recommend excluding at least the dates surrounding the event of interest.")
  }

  if(!identical(counts$date, arrange(counts,date)$date)) stop("counts must be ordered by date.")

  ## helper function
  fourier_trend <- function(x, k = 3){
    H <- lapply(1:k, function(k){
      cbind(sin(2*pi*k/365*x), cos(2*pi*k/365*x))
    })
    res <- do.call(cbind, H)
    colnames(res) <- paste(rep(c("sin", "cos"), k), rep(1:k, each = 2), sep="_")
    res
  }

  ## number of observations per year
  if(is.null(frequency)){
    frequency <- round(365 / (as.numeric(diff(range(counts$date)))/nrow(counts)))
    if(verbose) message("No frequency provided, determined to be ", frequency, " measurements per year.")
  }

  if(frequency < 365 & weekday.effect){
    warning("Modeling day effects is not recommended when frequency < 365. Consider setting weekday.effct = FALSE")
  }

  if(verbose) message("Overall death rate is ", signif(sum(counts$outcome, na.rm = TRUE)/sum(counts$population, na.rm = TRUE) * frequency * 1000, 3), ".")

  if(mean(counts$outcome, na.rm = TRUE) < 1)
    warning("Average counts per unit of time is below 1.")

  ## build design matrix
  # convert dates to time

  tt <- as.numeric(counts$date)
  index <- !(counts$date %in% exclude)

  # compute knots
  years <- (max(tt) - min(tt)) / 365
  nknots <- floor(years*trend.knots.per.year) + 1
  knots <- seq(min(tt), max(tt), length = nknots)
    
  # make trend basis (includes intercept)
  if(nknots > 2){
    
    knots <- knots[-c(1, length(knots))]
    x_t <- splines::ns(tt, knots = knots, intercept = TRUE)
    
  } else{
    
    knots <- c()
    x_t <- model.matrix(~tt)
  }
  
  # trend indices 
  i_t <- 1:ncol(x_t)

  #for harmonic model
  yd <- noleap_yday(counts$date)
  x_h <- fourier_trend(yd, k = harmonics)
  i_h <- ncol(x_t) + 1:ncol(x_h)

  ## build desing matrix
  if(weekday.effect){
    ## weekday effects
    w <- factor(lubridate::wday(counts$date))
    contrasts(w) <- contr.sum(length(levels(w)), contrasts = TRUE)

    x_w <- model.matrix(~w)[, -1] ## intercept already in spline
    i_w <- ncol(x_t) + ncol(x_h) + 1:ncol(x_w)
    x <- cbind(x_t, x_h, x_w)
  } else{
    x <- cbind(x_t, x_h)
  }

  y <- counts$outcome
  n <- counts$population

  ## fit model

  fit <- glm( y[index] ~ x[index,]-1, offset = log(n[index]), family = "quasipoisson")
  dispersion <- pmax(1, summary(fit)$dispersion)

  # prepare stuff to return
  expected <- exp(x %*% fit$coefficients) * n
  cova            <- summary(fit)$cov.unscaled * dispersion
  log_expected_se <- apply(X=x, MARGIN=1, function(xi){sqrt(t(xi) %*% cova %*% xi)})
  counts          <- mutate(counts, log_expected_se = as.numeric(log_expected_se))

  counts <- mutate(counts, expected = as.numeric(expected), log_expected_se = as.numeric(log_expected_se), excluded = !index)
  attr(counts, "dispersion") <- dispersion
  attr(counts, "trend.knots.per.year") <- trend.knots.per.year
  attr(counts, "knots") <- knots
  attr(counts, "harmonics") <- harmonics
  attr(counts, "frequency") <- frequency
  attr(counts, "weekday.effect") <- weekday.effect

  if(keep.components){

    seasonal <- data.frame(day = seq(1, 365, length = frequency),
                           s = exp(fourier_trend(seq(1, 365, length=frequency), k = harmonics)  %*% fit$coefficients[i_h]) -1)

    trend <- as.numeric(exp(x_t %*% fit$coefficients[i_t])  * frequency * 1000)

    if(weekday.effect){
      w <- factor(1:7)
      contrasts(w) <- contr.sum(length(levels(w)), contrasts = TRUE)
      weekday <- data.frame(weekday = 1:7,
                            effect = exp(model.matrix(~w)[, -1] %*% fit$coefficients[i_w])-1)
    } else{
      weekday <- NULL
    }

    res <- list(counts = counts,
                trend = trend,
                seasonal = seasonal,
                weekday = weekday)
  } else{
    res <- counts
  }

  attr(res, "class") <-  append(class(res), "compute_expected")
  attr(res, "keep.components") <- keep.components

  return(res)
}

