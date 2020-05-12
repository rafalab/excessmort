#' Fit excess count model
#' @export
#' @importFrom stats ARMAacf glm poly qnorm
excess_model <- function(counts,
                         event = NULL,
                         start = NULL,
                         end = NULL,
                         nknots = 12,
                         discontinuity = TRUE,
                         intervals = NULL,
                         model = c("quasipoisson", "poisson", "correlated"),
                         exclude = NULL,
                         knots.per.year = 1/5,
                         harmonics = 2,
                         frequency = NULL,
                         weekday.effect = TRUE,
                         control.dates = NULL,
                         max.control = 5000,
                         order.max = 14,
                         aic = TRUE,
                         max.iter = 15,
                         eps = 1e-8,
                         alpha = 0.05,
                         min.rate = 0.01,
                         verbose = TRUE){

  if("compute_expected" %in% class(counts)){
    if(attr(counts, "keep.components")){
      counts <- counts$counts
    }
  } else{
    if(verbose) message("Expected counts not provided, calculating them.")
    counts <-  compute_expected(counts,
                                exclude = exclude,
                                knots.per.year = knots.per.year,
                                harmonics = harmonics,
                                frequency = frequency,
                                weekday.effect = weekday.effect,
                                keep.components = FALSE)
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
    warning("Length of control longer than", max.control,". May result in long wait.")

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
    mu <- counts$expected[ind]
    obs <- counts$outcome[ind]
    pop <- counts$population[ind]

    ## compute residuals to fit ar model
    if(correlated.errors) y <- (obs - mu) / mu

    ## create the design matrix
    knots <- x[round(seq(1, n, length = nknots + 2))]
    knots <- knots[-c(1, length(knots))]
    if(!is.null(event)){
      event_index <- x[which(date == event)]
      i <- which.min(abs(knots - event_index))
      ##shift knots so that one of the internal knots falls on the event day
      knots <- knots + (event_index -  knots[i])
      X <- cbind(1, splines::ns(x, knots = knots))
      ## add parameters to account for discontinuity
      if(discontinuity){
        ind <- as.numeric(I(x >= event_index))
        X <- cbind(X, ind, poly((x - event_index)*ind, degree = 2))
      }
    } else{
      X <- cbind(1, splines::ns(x, knots = knots))
    }

    if(correlated.errors){
      fhat <- 0
      beta <- 0; beta0 <- 1
      count <- 0

      ## convinience function
      mysolve <- function(x) chol2inv(chol(x))

      ## parameters for covariance matrix
      if(length(arfit$ar) > 0 & arfit$sigma > 0){
        rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)
      }

      ## start iterations
      while(count < max.iter & sum((beta-beta0)^2) > eps){
        if(length(arfit$ar) > 0 & s > 0){
          Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) *
            outer(sqrt(s^2 + (1+fhat)/mu), sqrt(s^2 + (1+fhat)/mu))
          Sigma_inv <- mysolve(Sigma)
        } else{
          Sigma <- diag(s^2 + (1+fhat)/mu)
          Sigma_inv <- diag(1/(s^2 + (1+fhat)/mu))
        }
        ## fit spline using weighted least squares
        xwxi <- mysolve(t(X) %*% Sigma_inv %*% X)
        beta0 <- beta
        beta <- xwxi %*% t(X) %*% Sigma_inv %*% y
        count <- count + 1
        fhat <- pmax(as.vector(X %*% beta), min.rate - 1)
      }
      if(count >= max.iter) warning("No convergence after ", max.iter, " iterations.")

      se <- sqrt(apply(X, 1, function(x) matrix(x, nrow = 1) %*% xwxi %*% matrix(x, ncol = 1)))
      betacov <- xwxi

    } else{
      fit <- glm(obs ~ X-1, offset = log(mu), family = "poisson")
      tmp<- predict(fit, se = TRUE, type = "response")
      fhat <- tmp$fit/mu - 1
      se <- tmp$se * sqrt(dispersion) / mu
      Sigma <- diag(n) *  dispersion / mu
      betacov <- summary(fit)$cov.unscaled * dispersion
    }

    ## Compute regions for which estimate was outside usual range

    ind <- which(fhat - qnorm(1 - alpha/2) * se >= 0)
    if(length(ind) > 0){
      cluster <- cumsum(c(2, diff(ind)) > 1)
      indexes <- split(ind, cluster)
      excess <- lapply(indexes, function(ind){
        n <- length(ind)
        excess_stats(min(counts$date[ind]),
                     max(counts$date[ind]),
                     counts$outcome[ind],
                     counts$expected[ind],
                     Sigma[ind, ind, drop = FALSE],
                     counts$population[ind],
                     frequency,
                     fhat[ind],
                     X[ind,,drop=FALSE],
                     betacov)
      })
      detected_intervals <- do.call(rbind, excess)
    } else{
      detected_intervals <- data.frame(start = NA, end = NA, total = 0, se = NA, natural= NA)
    }
    ret <- list(date = date,
                observed = obs,
                expected = mu,
                fitted = fhat,
                population = pop,
                cov = Sigma,
                x = X,
                betacov = betacov,
                dispersion = dispersion,
                se = se,
                detected_intervals = detected_intervals)
    attr(ret, "frequency") <- frequency
    attr(ret, "model") <- match.arg(model)
    attr(ret, "type") <- "curve_fit"
    if(correlated.errors){
      ret$ar <-  arfit
    }
  } else{
    ret <- list()
  }

  ## If intervals provided compute excess deaths in each one
  if(!is.null(intervals)){
    res <- lapply(intervals, function(dates){
      ind <- which(counts$date %in% dates)
      date <- counts$date[ind]
      n <- length(ind)
      mu <- counts$expected[ind]

      if(correlated.errors){
        if(length(arfit$ar) > 0 & arfit$sigma > 0){
          rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)
        }
        if(length(arfit$ar) > 0 & arfit$sigma > 0){
          Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) *
            outer(sqrt(s^2 + 1/mu), sqrt(s^2 + 1/mu))
        } else{
          Sigma <- diag(s^2 + 1/mu)
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
    ret$excess <- do.call(rbind, res)
    attr(res, "type") <- append(attr(res, "type"), "excess")
  }

  return(ret)

}
