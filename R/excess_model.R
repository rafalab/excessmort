#' Fit excess count model
#' 
#' This function estimates the increase in the rate for a count time series relative to 
#' the rate for a typical year. Two options are available: 1 - model the rate increase as a 
#' smooth function and estimate this function or 2 - estimate the total excess in intervals. 
#' For 1m an `event` date can be provided and a discontinuity included in the model.
#' You can do either 1 or 2 or both. 
#' 
#' Three versions of the model are available: 1 - Assume counts are Poisson distributed, 
#' 2 - assume counts are overdispersed Poisson, or 3 - assume a mixed model with 
#' correlated errors. The second is the default and recommended for weekly count data. For daily counts we often find 
#' evidence of correlation and recommend the third along with setting `weekday.effect = TRUE`.
#' 
#' If the `counts` object includes a `expected` column produced by `compute_expected` these are used
#' as the expected counts. If not, then these are computed.
#' 
#' 
#' @param counts A data frame with date, count and population columns.
#' @param start First day of interval to which model will be fit
#' @param end Last day of interval to which model will be fit
#' @param knots.per.year Number of knots per year used for the fitted smooth function
#' @param event If modeling a discontinuity is desired, this is the day in which it happens
#' @param intervals Instead of `start` and `end` a list of time intervals can be provided and excess is computed in each one
#' @param discontinuity Logical that determines if discontinuity is allowed at `event`
#' @param model Which version of the model to fit
#' @param exclude Dates to exclude when computing expected counts
#' @param trend.knots.per.year Number of knots per year used by `compute_expected` to estimate the trend for the expected counts
#' @param harmonics Number of harmonics used by `compute_expected` to estimate seasonal trend
#' @param frequency Number of observations per year. If not provided an attempt is made to calculate it
#' @param weekday.effect Logical that determins if a day of the week effects is included in the model. Should be `FALSE` for weekly or monthly data.
#' @param control.dates When `model` is set to `correlated`, these dates are used to estimate the covariance matrix. The larger this is the slower the function runs.
#' @param max.control If the length of `control.dates` is larger than `max.control` the function stops.
#' @param order.max Larges order for the Autoregressive process used to model the covariance structure
#' @param aic A logical that determines if the AIC criterion is used to selected the order of the AR process
#' @param maxit Maxium number of iterations for the IRLS algorithm used when `model` is `correlated`
#' @param epsilon Difference in deviance requried to declare covergenace of IRLS
#' @param alpha Percentile used to define what is outside the normal range
#' @param min.rate The estimated expected rate is not permited to go below this value
#' @param verbose Logical that determines if messages are displayed
#' 
#' @return If only `intervals` are provided a data frame with excess estimates described below for `excess`. 
#' if `start` and `end` are provided the a list with the following components are included:
#' \describe{
#' \item{date}{The dates for which the estimate was computed}
#' \item{observed}{The observed counts}
#' \item{expected}{The expected counts}
#' \item{fitted}{The fitted curve for excess counts}
#' \item{se}{The point-wise standard error for the fitted curve}
#' \item{population}{The population size}
#' \item{sd}{The standard deviation for observed counts on a typical year}
#' \item{cov}{The estimated covariance matrix for the observed counts}
#' \item{x}{The design matrix used for the fit}
#' \item{betacov}{The covariance matrix for the estimated coefficients}
#' \item{dispersion}{The estimated overdispersion parameter}
#' \item{detected_intervals}{Time intervals for which the 1 - `alpha` confidence interval does not include 0}
#' \item{ar}{The estimated coefficients for the autoregressive process}
#' \item{excess}{A data frame with information for the time intervals provided in `itervals`. This includes start, end, observed death rate (per 1,000 per year), expected death rate, standard deviation for the death rate, observed counts, expected counts, excess counts, standard deviation}
#' }
#' 
#' @examples
#' data(cdc_state_counts)
#' counts <- cdc_state_counts[cdc_state_counts$state == "Massachusetts", ]
#' exclude_dates <- c(seq(as.Date("2017-12-16"), as.Date("2018-01-16"), by = "day"),
#' seq(as.Date("2020-01-01"), max(cdc_state_counts$date), by = "day"))
#' f <- excess_model(counts, 
#' exclude = exclude_dates,
#' start = min(counts$date),
#' end = max(counts$date),
#' knots.per.year = 12)

#' data(new_jersey_counts)
#' exclude_dates <- as.Date("2012-10-29") + 0:180
#' control_dates <- seq(min(new_jersey_counts$date), min(exclude_dates) - 1, by = "day")
#' f <- excess_model(new_jersey_counts,
#' start = as.Date("2012-09-01"), 
#' end = as.Date("2013-09-01"), 
#' exclude = exclude_dates,
#' model = "correlated",
#' weekday.effect = TRUE,
#' control.dates = control_dates)
#' 
#' @export
#' @importFrom stats ARMAacf glm poly qnorm fitted.values
excess_model <- function(counts,
                         start = NULL,
                         end = NULL,
                         knots.per.year = 12,
                         event = NULL,
                         intervals = NULL,
                         discontinuity = TRUE,
                         model = c("quasipoisson", "poisson", "correlated"),
                         exclude = NULL,
                         trend.knots.per.year = 1/7,
                         harmonics = 2,
                         frequency = NULL,
                         weekday.effect = FALSE,
                         control.dates = NULL,
                         max.control = 5000,
                         order.max = 14,
                         aic = TRUE,
                         maxit = 25,
                         epsilon = 1e-8,
                         alpha = 0.05,
                         min.rate = 0.0001,
                         verbose = TRUE){

  if("compute_expected" %in% class(counts)){
    if(attr(counts, "keep.components")){
      counts <- counts$counts
    }
  } else{
    if(verbose) message("Computing expected counts.")
    counts <-  compute_expected(counts,
                                exclude = exclude,
                                trend.knots.per.year = trend.knots.per.year,
                                harmonics = harmonics,
                                frequency = frequency,
                                weekday.effect = weekday.effect,
                                keep.components = FALSE,
                                verbose = verbose)
  }

  correlated.errors <- match.arg(model) == "correlated"

  ## number of observations per year
  frequency <- attr(counts, "frequency")
  dispersion <- attr(counts, "dispersion")

  ## checks
  if(any(counts$excluded)) exclude <- counts$date[counts$excluded] else exclude <- NULL

  if(correlated.errors & is.null(control.dates)){
    if(!is.null(exclude)){
      warning("No control region suplied, which is not recommended when correlated.errors = TRUE. Using data up to first excluded point.")
      control.dates <-seq(min(counts$date), min(exclude, na.rm = TRUE) - 1, by = "day")
      } else{
      warning("No control region suplied, which is not recommended when correlated.errors = TRUE. Using all the data")
      control.dates <- counts$date
    }
  }

  if(length(control.dates) > max.control & correlated.errors)
   stop("Length of control longer than", max.control)

  if((is.null(start) & !is.null(end)) | (!is.null(start) & is.null(end)))
    stop("You must provide both start and end, not just one.")

  if(is.null(start) & is.null(end) & is.null(intervals))
    stop("You must provide start and end or intervals.")


  ## check to see if counts per unit of time are high enough for model to work
  if(mean(counts$outcome, na.rm = TRUE) < 1 & correlated.errors)
    warning("Low counts per unit of time, consider fitting model with no correlation.")

  ## compute_expected always uses quasipoisson
  if(match.arg(model) == "poisson") dispersion <- 1

  ## Use control days to compute the autocorrelation function
  if(correlated.errors){
    arfit <- fit_ar(counts, control.dates, order.max = order.max, aic = aic)
    s <- arfit$sigma
    if(verbose) message("Order selected for AR model is ",length(arfit$ar),
                        ". Estimated residual standard error is ", signif(arfit$sigma, 2), ".")
  }

  ## Fit start and end provided, fit the curve
  if(!is.null(start) & !is.null(end)){

    if(!is.null(event)){
      if(event >= end | event <= start)
        stop("event must be between start and end.")
    }

    ## now fit the GLS to the relevant subset of data
    include_dates = seq(start, end, by = "day")

    if(frequency == 365)
      include_dates <- include_dates[!(lubridate::month(include_dates) == 2 & lubridate::day(include_dates) == 29)]

    ind <- which(counts$date %in% include_dates)
    date <- counts$date[ind]
    n <- length(ind)
    x <- 0:(n-1) / frequency * 365
    mu <- pmax(min.rate, counts$expected[ind])
    if(any(mu == min.rate)) warning("Minimum expected rate reached and was set at ", min.rate, ".")
    min.fhat <- min.rate / mu - 1
    obs <- counts$outcome[ind]
    pop <- counts$population[ind]

    ## compute residuals to fit ar model
    if(correlated.errors) y <- (obs - mu) / mu

    ## create the design matrix
    nknots <- round(knots.per.year * as.numeric(max(date) - min(date)) / 365)
    knots <- x[round(seq(1, n, length = nknots + 2))]
    knots <- knots[-c(1, length(knots))]
    if(!is.null(event)){
      event_index <- x[which.min(abs(as.numeric(date - event)))]
      i <- which.min(abs(knots - event_index))
      ##shift knots so that one of the internal knots falls on the event day
      knots <- knots + (event_index -  knots[i])
      X <- cbind(1, splines::ns(x, knots = knots))
      ## add parameters to account for discontinuity
      if(discontinuity){
        after_ind <- as.numeric(I(x >= event_index))
        X <- cbind(X, after_ind, poly((x - event_index)*after_ind, degree = 2))
      }
    } else{
      X <- cbind(1, splines::ns(x, knots = knots))
    }

    if(correlated.errors){
      fhat <- 0
      beta <- 0;
      dev<- 2*sum(ifelse(obs == 0, 0, obs * log(obs / mu)) - (obs - mu))
      ## convinience function
      mysolve <- function(x) chol2inv(chol(x))

      ## parameters for covariance matrix
      if(length(arfit$ar) > 0 & arfit$sigma > 0){
        rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)
      }

      ## start iterations
      count <- 0
      flag <- TRUE
      log_mu_vari <- counts$log_expected_se[ind]^2
      while(count < maxit & flag){
        if(length(arfit$ar) > 0 & s > 0){
          Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) * 
            outer(sqrt((1 + fhat)^2 * s^2 + (1 + fhat)/mu + (1 + fhat)^2 * log_mu_vari), sqrt((1 + fhat)^2 * s^2 + (1 + fhat)/mu + (1 + fhat)^2 * log_mu_vari))
          Sigma_inv <- mysolve(Sigma)
        } else{
          Sigma <- diag((1 + fhat)^2 * s^2 + (1 + fhat)/mu + (1 + fhat)^2 * log_mu_vari)
          Sigma_inv <- diag(1/((1 + fhat)^2 * s^2 + (1 +  fhat)/mu + (1 + fhat)^2 * log_mu_vari))
        }
        ## fit spline using weighted least squares
        xwxi <- mysolve(t(X) %*% Sigma_inv %*% X)
        beta <- xwxi %*% t(X) %*% Sigma_inv %*% y
        count <- count + 1
        fhat <- pmax(as.vector(X %*% beta), min.fhat)
        devold <- dev
        dev <- 2*sum(ifelse(obs == 0, 0, obs * log(obs / (mu*(1 + fhat)))) - (obs - mu*(1 + fhat)))
        flag <- abs(dev - devold)/(0.1 + abs(dev)) >= epsilon
      }
      if(count > maxit) warning("No convergence after ", maxit, " imterations.")

      se <- sqrt(apply(X, 1, function(x) matrix(x, nrow = 1) %*% xwxi %*% matrix(x, ncol = 1)))
      betacov <- xwxi

    } else{
      fit <- glm(obs ~ X-1, offset = log(mu), family = "poisson")
      tmp <- predict(fit, se = TRUE, type = "response")
      fhat <- pmax(tmp$fit/mu - 1, min.fhat)
      lambda <- fitted.values(fit)
      cova <- summary(fit)$cov.unscaled * summary(fit)$dispersion
      lambda_vari <- lambda^2 * diag(X %*% cova %*% t(X))
      mu_vari <- counts$expected[ind]^2 * counts$log_expected_se[ind]^2
      se <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
      Sigma <- diag(n) *  summary(fit)$dispersion / mu
      betacov <- summary(fit)$cov.unscaled * summary(fit)$dispersion
    }

    ## Warning if minimum reached
    if(any(fhat == min.fhat)) warning("Minimum rate reached and was set so that estimated rate is ", min.rate, ".")
    
    ## Compute regions for which estimate was outside usual range

    excess_ind <- which(fhat - qnorm(1 - alpha/2) * se >= 0)
    if(length(excess_ind) > 0){
      cluster <- cumsum(c(2, diff(excess_ind)) > 1)
      indexes <- split(excess_ind, cluster)
      excess <- lapply(indexes, function(i){
        n <- length(i)
        excess_stats(min(counts$date[ind[i]]),
                     max(counts$date[ind[i]]),
                     counts$outcome[ind[i]],
                     counts$expected[ind[i]],
                     Sigma[i, i, drop = FALSE],
                     counts$population[ind[i]],
                     frequency,
                     fhat[i],
                     X[i,,drop=FALSE],
                     betacov)
      })
      detected_intervals <- do.call(rbind, excess)
    } else{
      detected_intervals <- data.frame(start = NA, end = NA, total = 0, se = NA, natural= NA)
    }
    ret <- list(date = date,
                observed = obs,
                expected = mu,
                log_expected_se = counts$log_expected_se[ind],
                fitted = fhat,
                se = se,
                population = pop,
                sd = sqrt(diag(Sigma)),
                cov = Sigma,
                x = X,
                betacov = betacov,
                dispersion = dispersion,
                detected_intervals = detected_intervals)

    attr(ret, "frequency") <- frequency
    attr(ret, "model") <- match.arg(model)
    attr(ret, "type") <- "curve_fit"

    if(correlated.errors){
      ret$ar <-  arfit
    }

  }

  ## If intervals provided compute excess deaths
  ## The uncertainty is calculated under the null
  if(!is.null(intervals)){
    res <- lapply(intervals, function(dates){
      ind         <- which(counts$date %in% dates)
      date        <- counts$date[ind]
      n           <- length(ind)
      mu          <- counts$expected[ind]
      log_mu_vari <- counts$log_expected_se[ind]^2

      if(correlated.errors){
        if(length(arfit$ar) > 0 & arfit$sigma > 0){
          rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)
        }
        if(length(arfit$ar) > 0 & arfit$sigma > 0){
          Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) *
            outer(sqrt(s^2 + 1/mu + log_mu_vari), sqrt(s^2 + 1/mu + log_mu_vari))
        } else{
          Sigma <- diag(s^2 + 1/mu + log_mu_vari)
        }
      } else{
        Sigma <- diag(n) *  dispersion / mu
      }

      excess_stats(min(dates),
                   max(dates),
                   counts$outcome[ind],
                   mu,
                   Sigma,
                   counts$population[ind],
                   frequency)
    })
    res <- do.call(rbind, res)
    if(!exists("ret")){
      ret <- res
      attr(ret, "type") <- "excess"
    }  else{
      ret$excess <- res
      attr(ret, "type") <- append(attr(ret, "type"), "excess")
    }
  }
  return(ret)
}

