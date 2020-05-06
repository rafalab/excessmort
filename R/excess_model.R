#' Fit excess count model
#' @export
#' @importFrom stats ARMAacf glm poly qnorm
excess_model <- function(counts,
                         event = NULL,
                         start = NULL,
                         end = NULL,
                         expected = NULL,
                         exclude = NULL,
                         trend.nknots = 1/5,
                         harmonics = 2,
                         family = "quasipoisson",
                         day.effect = TRUE,
                         nknots = 4,
                         discontinuity = TRUE,
                         max.iter = 15,
                         eps = 1e-8,
                         alpha = 0.05,
                         min.rate = 0.01){ ##smallest possible death rate

  ## order by dates
  counts <- counts %>% arrange(date)

  ## checks
  if(is.null(event) & (is.null(start) | is.null(end))){
    warning("No event or start/end provide, using all data")
    start <- min(counts$date)
    end <- max(counts$date)
  }
  if(is.null(start)) start = event - lubridate::years(1)
  if(is.null(end)) end = event + lubridate::years(1)

  if(mean(counts$outcome, na.rm = TRUE) < 1)
    warning("Average counts per unit of time below 1.")

  if(is.null(exclude)) warning("No dates excluded.")
  ## number of observations per year
  TT <- round(365 / (as.numeric(diff(range(counts$date)))/nrow(counts)))

  ## if expected not supplied compute it
  if(is.null(expected)){
    expected <- compute_expected(counts, exclude = exclude,
                                 trend.nknots = trend.nknots,
                                 harmonics = harmonics,
                                 family = family,
                                 day.effect = day.effect)
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

  fit <- glm(obs ~ X-1, offset = log(mu), family = "poisson")
  tmp<- predict(fit, se = TRUE)
  fhat <- tmp$fit - log(mu)
  se <- tmp$se * sqrt(expected$dispersion)

  ind <- which(fhat - qnorm(1 - alpha/2)*se >= 0)
  if(length(ind) > 0){
    cluster <- cumsum(c(2, diff(ind)) > 1)
    indexes <- split(ind, cluster)
    excess <- lapply(indexes, function(ind){
      excess <- mu[ind] %*% fhat[ind]

      excess_se <- sqrt(matrix(mu[ind], nrow = 1) %*%
                          X[ind,,drop=FALSE] %*%
                          summary(fit)$cov.unscaled%*%
                          t(X[ind,,drop=FALSE]) %*%
                          matrix(mu[ind], ncol = 1))
      excess_se <- excess_se * sqrt(expected$dispersion)

      natural_se <- mu[ind] * sqrt(expected$dispersion)
      data.frame(start = date[ind[1]], end = date[ind[length(ind)]], total = excess,  se = excess_se, natural = natural_se)
    })
    excess <- do.call(rbind, excess)
  } else excess <- data.frame(start = NA, end = NA, total = 0, se = NA, natural= NA)
  return(list(date = date,
              observed = obs,
              expected = mu,
              fitted = fhat,
              population = pop,
              x = X,
              betacov = summary(fit)$cov.unscaled*sqrt(expected$dispersion),
              se = se,
              excess = excess,
              plugins = expected
              ))
}
