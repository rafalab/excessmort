excess_model <- function(counts,
                         event = NULL,
                         start = NULL,
                         end = NULL,
                         control.dates = NULL,
                         expected = NULL,
                         exclude = NULL,
                         trend.nknots = 1/5,
                         harmonics = 2,
                         family = "poisson",
                         max.control = 5000,
                         nknots = 4,
                         discontinuity = TRUE,
                         order.max = 14,
                         aic = TRUE,
                         max.iter = 15,
                         eps = 1e-8,
                         alpha = 0.05,
                         min.rate = 0.01){ ##smallest possible death rate

  if(is.null(event) & (is.null(start) | is.null(end))) stop("Need to provide event or start/end")

  if(is.null(start)) start = event - years(1)
  if(is.null(end)) end = event + years(1)

  if(mean(counts$outcome, na.rm = TRUE) < 1)
    warning("Average counts per unit of time below 1.")

  if(is.null(control.dates)){
    warning("No control region suplied, using all data. This is not recommended.")
    control.dates <- expected$date
  }

  if(length(control.dates) > max.control) warning("Length of control longer than", max.control,". May result in long wait.")

  if(is.null(exclude)) warning("No dates excluded.")
  ## number of observations per year
  TT <- length(unique(noleap_yday(counts$date)))

  ## if expected not supplied compute it
  if(is.null(expected)){
    expected <- compute_expected(counts, exclude = exclude,
                                 trend.nknots = trend.nknots,
                                 harmonics = harmonics,
                                 family = family)
  }

  ## Use control days to compute the autocorrelation function
  arfit <- fit_ar(expected, control.dates, order.max = order.max, aic = aic)

  ## now fit the GLS to the relevant subset of data
  include_dates = seq(start, end, by = "day")
  include_dates <- include_dates[!(month(include_dates) == 2 & day(include_dates) == 29)]
  ind <- which(expected$date %in% include_dates)
  date <- expected$date[ind]
  n <- length(ind)
  x <- 1:n
  y <- expected$resid[ind]
  mu <- expected$expected[ind]
  obs <- counts$outcome[ind]

  ## create the design matrix
  df <- round(n / TT * nknots)
  knots <- round(seq(1, n, length = df))
  knots <- knots[-c(1, length(knots))]
  if(!is.null(event)){
    event_index <- x[which(date == event)]
    i <- which.min(abs(knots - event_index))
    ##shift knots so that one of the internal knots falls on the event day
    knots <- knots + (event_index -  knots[i])
    X <- cbind(1, ns(x, knots = knots))
    ## add parameters to account for discontinuity
    if(discontinuity){
      ind <- as.numeric(I(x >= event_index))
      X <- cbind(X, ind, poly((x - event_index)*ind, degree = 2))
    }
  } else{
    X <- cbind(1, ns(x, knots = knots))
  }

  ## Fit using interative algorithm
  ## starting values
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
    if(length(arfit$ar) > 0 & arfit$sigma > 0){
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

  ind <- which(fhat - qnorm(1 - alpha/2)*se >= 0)
  if(length(ind) > 0){
    cluster <- cumsum(c(2, diff(ind)) > 1)
    indexes <- split(ind, cluster)
    excess <- map_df(indexes, function(ind){
      excess <- mu[ind] %*% fhat[ind]

      excess_se <- sqrt(matrix(mu[ind], nrow = 1) %*%
                          X[ind,,drop=FALSE] %*%
                          xwxi %*%
                          t(X[ind,,drop=FALSE]) %*%
                          matrix(mu[ind], ncol = 1))

      natural_se <- sqrt(matrix(mu[ind], nrow = 1) %*% Sigma[ind, ind] %*% matrix(mu[ind], ncol= 1))
      list(start = date[ind[1]], end = date[ind[length(ind)]], total = excess,  se = excess_se, natural = natural_se)
    })
  } else excess <- data.frame(start = NA, end = NA, total = 0, se = NA, natural= NA)
  return(list(date = date,
              observed = obs,
              expected = mu,
              resid = y,
              trend = expected$trend,
              seasonal = expected$seasonal,
              weekday = expected$weekday,
              x = X,
              betacov = xwxi,
              fitted = fhat,
              se = se,
              excess = excess,
              ar = arfit,
              cov = Sigma))
}
