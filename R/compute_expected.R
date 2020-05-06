#' Compute expected counts for each day
#' @export
#' @importFrom stats glm contr.sum model.matrix contrasts<-
#'
compute_expected <- function(counts, exclude = NULL,
                             trend.nknots = 1/5,
                             harmonics = 2,
                             family = "quasipoisson",
                             day.effect = TRUE){

  ## helper function
  fourier_trend <- function(x, k = 3){
    H <- lapply(1:k, function(k){
      cbind(sin(2*pi*k/365*x), cos(2*pi*k/365*x))
    })
    res <- do.call(cbind, H)
    colnames(res) <- paste(rep(c("sin", "cos"), k), rep(1:k, each = 2), sep="_")
    res
  }
  ## order by dates
  counts <- counts <- counts %>% arrange(date)

  ## number of observations per year
  TT <- round(365 / (as.numeric(diff(range(counts$date)))/nrow(counts)))

  ## build design matrix
  # compute dfs
  dfs <- round(length(unique(lubridate::year(counts$date[!counts$date %in% exclude])))*trend.nknots)

  # make trend basis (includes intercept)
  x_t <- splines::ns(as.numeric(counts$date), df = dfs + 1, intercept = TRUE)
  i_t <- 1:ncol(x_t)

  #for harmonic model
  yd <- noleap_yday(counts$date)
  x_h <- fourier_trend(yd, k = harmonics)
  i_h <- ncol(x_t) + 1:ncol(x_h)

  ## build desing matrix
  if(day.effect){
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
  index <- which(!counts$date %in% exclude)

  fit <- glm( y[index] ~ x[index,]-1, offset = log(n[index]), family = family)
  dispersion <- pmax(1, summary(fit)$dispersion)
  # prepare stuff to return
  expected <- exp(x %*% fit$coefficients) * n

  seasonal <- data.frame(day = seq(0, 364, length=TT),
                         s = exp(fourier_trend(seq(0, 364, length=TT), k = harmonics)  %*% fit$coefficients[i_h]) -1)

  trend <- exp(x_t %*% fit$coefficients[i_t])  * TT * 1000

  if(day.effect){
    w <- factor(1:7)
    contrasts(w) <- contr.sum(length(levels(w)), contrasts = TRUE)
    weekday <- data.frame(weekday = 1:7,
                          effect = exp(model.matrix(~w)[, -1] %*% fit$coefficients[i_w])-1)
  } else{
    weekday <- NULL
  }

  ## add expected counts to data table
  return(list(date = counts$date,
              observed = counts$outcome,
              expected = expected,
              dispersion = dispersion,
              trend = trend,
              seasonal = seasonal,
              weekday = weekday))
}
