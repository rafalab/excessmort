#' Helper funtion to compute excess deaths statistics for a
#' @importFrom stats ARMAacf
#'
#'
excess_stats <- function(start, end, obs, mu, cov, pop, frequency,
                         fhat = NULL, X = NULL, betacov = NULL){

  mu <- matrix(mu, nrow = 1)

  observed <- sum(obs)

  expected <- sum(mu)

  excess <- observed - expected

  ##  cov is the variance-covariance of the percent increase                                                                                                                                                                                                              (obs - exp) / exp
  excess_se <- sqrt(mu %*% cov %*% t(mu))

  res <- data.frame(
    start = start,
    end = end,
    obs_death_rate = observed / sum(pop) * frequency * 1000,
    exp_death_rate = expected / sum(pop) * frequency * 1000,
    se_death_rate = excess_se / sum(pop) * frequency * 1000,
    observed = observed,
    expected = expected,
    excess = excess,
    excess_se = excess_se)

  if(!is.null(fhat)){
    res$fitted <- mu %*% fhat
    res$fitted_se <- sqrt(mu %*% X %*% betacov %*% t(X) %*% t(mu))
  }

  return(res)
}

