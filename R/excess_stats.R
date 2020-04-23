excess_stats <- function(counts,
                         intervals,
                         control.dates,
                         expected = NULL,
                         exclude = NULL,
                         trend.nknots = 1/5,
                         harmonics = 2,
                         family = "poisson",
                         max.control = 5000,
                         order.max = 30,
                         aic = TRUE){

  if(is.null(control.dates)){
    warning("No control region suplied, using all data. This is not recommended.")
    control.dates <- expected$date
  }

  if(length(control.dates) > max.control) warning("Length of control longer than", max.control,". May result in long wait.")

  if(is.null(exclude)) warning("No dates excluded.")

  ## if expected not supplied compute it
  if(is.null(expected)){
    expected <- compute_expected(counts, exclude = exclude,
                                 trend.nknots = trend.nknots,
                                 harmonics = harmonics,
                                 family = family)
  }


  ## Use control days to compute the autocorrelation function
  arfit <- fit_ar(expected, control.dates, order.max = order.max, aic = aic)
  s <- arfit$sigma

  date <- counts$date

  map_df(intervals, function(dates){
    ind <- which(date %in% dates)
    n <- length(ind)
    y <- expected$resid[ind]
    mu <- expected$expected[ind]
    obs <- counts$outcome[ind]
    pop <- counts$population[ind]
    date <- date[ind]

    if(length(arfit$ar) > 0 & arfit$sigma > 0){
      rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)
    }
    if(length(arfit$ar) > 0 & arfit$sigma > 0){
      Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) *
        outer(sqrt(s^2 + 1/mu), sqrt(s^2 + 1/mu))
    } else{
      Sigma <- diag(s^2 + 1/mu)
    }

    total <- sum(obs)
    expected <- sum(mu)
    se <- sqrt(matrix(mu, nrow = 1) %*%  Sigma %*%  matrix(mu, ncol = 1))
    return(list(start = date[1], end = date[length(ind)],
                observed = total, expected = expected, se = se,
                observed_death_rate = total / sum(pop) * 365 * 1000,
                expected_death_rate = expected / sum(pop) * 365 * 1000,
                se_death_rate = se / sum(pop) * 365 * 1000))
  })
}
