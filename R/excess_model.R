#' Fit excess count model
#' @export
#' @importFrom stats ARMAacf glm poly qnorm
excess_model <- function(counts,
                         event = NULL,
                         start = NULL,
                         end = NULL,
                         nknots = 12,
                         discontinuity = TRUE,
                         expected = NULL,
                         exclude = NULL,
                         trend.nknots = 1/5,
                         harmonics = 2,
                         day.effect = TRUE,
                         correlated.errors = FALSE,
                         control.dates = NULL,
                         max.control = 5000,
                         order.max = 14,
                         aic = TRUE,
                         max.iter = 15,
                         eps = 1e-8,
                         alpha = 0.05,
                         min.rate = 0.01,
                         verbose = TRUE){

  ## checks
  if(!identical(counts, arrange(counts,date))) stop("counts must be ordered by date.")

  if(any(table(counts$date))>1) stop("Each date can appeared at most once.")

  if(any(is.na(counts$date)) | any(is.na(counts$outcome)) | any(is.na(population)))
    stop("No NAs permited in date, outcome, or population columns.")

  ## number of observations per year
  TT <- round(365 / (as.numeric(diff(range(counts$date)))/nrow(counts)))
  if(verbose) message("Detected ",TT," measurements per year.")

  if(is.null(exclude)){
    warning("No dates excluded. We recommend excluding at least the dates surrounding the event of interest.")
  }

  if(correlated.errors & is.null(control.dates) & !is.null(exclude)){
    warning("No control region suplied, which is not recommended when correlated.errors = TRUE. Using data up to first excluded point.")
    control.dates <-seq(min(counts$date), min(exclude, na.rm = TRUE) - 1)
  }

  if(correlated.errors & is.null(control.dates) & is.null(exclude)){
    warning("No control region suplied, which is not recommended when correlated.errors = TRUE. Using all the data")
    control.dates <- counts$date
  }

  if(length(control.dates) > max.control & correlated.errors)
    warning("Length of control longer than", max.control,". May result in long wait.")

  ## check if event, start, or end are provided
  if(is.null(event) & (is.null(start) | is.null(end))){
    warning("No event or start/end provide, using all data")
    start <- min(counts$date)
    end <- max(counts$date)
  }
  ## event is provided with start and end, add a 1 year around event
  if(is.null(start)) start = event - lubridate::years(1)
  if(is.null(end)) end = event + lubridate::years(1)

  if(verbose) message("Overall death rate is ", signif(sum(counts$outcome, na.rm = TRUE)/sum(counts$population, na.rm = TRUE) * TT * 1000, 3), ".")

  ## check to see if counts per unit of time are high enough for model to work
  if(mean(counts$outcome, na.rm = TRUE) < 1 & correlated.errors)
    warning("Average counts per unit of time below 1. Consider fitting model with no correlation.")
  if(mean(counts$outcome, na.rm = TRUE) < 1 & !correlated.errors)
    warning("Average counts per unit of time below 1.")

  ## if expected not supplied compute it
  if(is.null(expected)){
    if(verbose) message("Expected counts not provided. Calculating them.")
    expected <- compute_expected(counts,
                                 exclude = exclude,
                                 trend.nknots = trend.nknots,
                                 harmonics = harmonics,
                                 day.effect = day.effect)
  }

  ## Use control days to compute the autocorrelation function
  if(correlated.errors){
    arfit <- fit_ar(expected, control.dates, order.max = order.max, aic = aic)
    if(verbose) message("Order selected fpr AR model is ",length(arfit$ar),
                        ". Estimated residual standard error is ", signif(arfit$sigma,2), ".")
  }

  ## now fit the GLS to the relevant subset of data
  include_dates = seq(start, end, by = "day")
  include_dates <- include_dates[!(lubridate::month(include_dates) == 2 & lubridate::day(include_dates) == 29)]
  ind <- which(expected$date %in% include_dates)
  date <- expected$date[ind]
  n <- length(ind)
  x <- 0:(n-1) / TT * 365
  mu <- expected$expected[ind]
  obs <- counts$outcome[ind]
  pop <- counts$population[ind]

  ## compute residuals to fit ar model
  if(correlated.errors) y <- (obs - mu) / mu

  ## create the design matrix
  df <- round(n / TT * nknots)
  knots <- x[round(seq(1, n, length = df))]
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
    s <- arfit$sigma
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

  } else{
    fit <- glm(obs ~ X-1, offset = log(mu), family = "poisson")
    tmp<- predict(fit, se = TRUE, type = "response")
    fhat <- tmp$fit/mu - 1
    se <- tmp$se * sqrt(expected$dispersion) / mu
  }

  ind <- which(fhat - qnorm(1 - alpha/2) * se >= 0)
  if(length(ind) > 0){
    cluster <- cumsum(c(2, diff(ind)) > 1)
    indexes <- split(ind, cluster)
    excess <- lapply(indexes, function(ind){
      excess <- mu[ind] %*% fhat[ind]

      if(correlated.errors){
        excess_se <- sqrt(matrix(mu[ind], nrow = 1) %*%
                            X[ind,,drop=FALSE] %*%
                            xwxi %*%
                            t(X[ind,,drop=FALSE]) %*%
                            matrix(mu[ind], ncol = 1))

        natural_se <- sqrt(matrix(mu[ind], nrow = 1) %*% Sigma[ind, ind] %*% matrix(mu[ind], ncol= 1))
      } else{
        excess_se <- sqrt(matrix(mu[ind], nrow = 1) %*%
                            X[ind,,drop=FALSE] %*%
                            summary(fit)$cov.unscaled %*%
                            t(X[ind,,drop=FALSE]) %*%
                            matrix(mu[ind], ncol = 1))
        excess_se <- excess_se * sqrt(expected$dispersion)
        natural_se <- sqrt(sum(mu[ind] * expected$dispersion))
      }

      data.frame(start = date[ind[1]], end = date[ind[length(ind)]], total = excess,  se = excess_se, natural = natural_se)
    })
    excess <- do.call(rbind, excess)
  } else excess <- data.frame(start = NA, end = NA, total = 0, se = NA, natural= NA)

  ## return results
  return(list(date = date,
              observed = obs,
              expected = mu,
              fitted = fhat,
              population = pop,
              x = X,
              betacov = if(correlated.errors) xwxi else summary(fit)$cov.unscaled*sqrt(expected$dispersion),
              se = se,
              frequency = TT,
              excess = excess,
              plugins = expected,
              correlated.errors = correlated.errors,
              cov = if(correlated.errors) Sigma else  NULL,
              ar = if(correlated.errors) arfit else NULL
  ))
}
