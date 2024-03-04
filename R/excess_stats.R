#' Excess counts in an interval
#' 
#' Helper function to compute excess deaths statistics for a
#' 
#' @param start First day of interval
#' @param end Last day of interval
#' @param obs Observed counts
#' @param mu Expected counts
#' @param cov Covariance matrix for percent change
#' @param pop Population size
#' @param frequency Observations per year
#' @param fhat Estimated percent increase
#' @param X Design matrix used to estimate fhat
#' @param betacov Covariance matrix for parameter estimates used to estimate fhat
#' 
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
  sd <- sqrt(mu %*% cov %*% t(mu))

  res <- data.frame(
    start = start,
    end = end,
    obs_death_rate = observed / sum(pop) * frequency * 1000,
    exp_death_rate = expected / sum(pop) * frequency * 1000,
    sd_death_rate = sd / sum(pop) * frequency * 1000,
    observed = observed,
    expected = expected,
    excess = excess,
    sd = sd)

  if(!is.null(fhat)){
    res$fitted <- mu %*% fhat
    res$se <- sqrt(mu %*% X %*% betacov %*% t(X) %*% t(mu))
  }

  return(res)
}

