#' Compute excess deaths statistics for a list of intervals
#' @export
#' @importFrom stats ARMAacf
excess_stats <- function(counts,
                         intervals,
                         expected = NULL,
                         exclude = NULL,
                         trend.nknots = 1/5,
                         harmonics = 2,
                         family = "quasipoisson",
                         day.effect = TRUE,
                         max.control = 5000){


  if(is.null(exclude)) warning("No dates excluded.")

  ## if expected not supplied compute it
  if(is.null(expected)){
    expected <- compute_expected(counts, exclude = exclude,
                                 trend.nknots = trend.nknots,
                                 harmonics = harmonics,
                                 family = family,
                                 day.effect = day.effect)
  }


  ## Use control days to compute the autocorrelation function
  date <- counts$date

  res <- lapply(intervals, function(dates){
    ind <- which(date %in% dates)
    n <- length(ind)
    mu <- expected$expected[ind]
    obs <- counts$outcome[ind]
    pop <- counts$population[ind]
    date <- date[ind]

    s2 <- mu*expected$dispersion

    total <- sum(obs)
    expected <- sum(mu)
    se <- sqrt(sum(s2))
    TT <- round(365 / (as.numeric(diff(range(counts$date)))/nrow(counts)))

    return(data.frame(start = date[1], end = date[length(ind)],
                      observed = total, expected = expected, se = se,
                      observed_death_rate = total / sum(pop) * TT * 1000,
                      expected_death_rate = expected / sum(pop) * TT * 1000,
                      se_death_rate = se / sum(pop) * TT * 1000))
  })

  res <- do.call(rbind, res)

  return(res)
}
