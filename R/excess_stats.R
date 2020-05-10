#' Compute excess deaths statistics for a list of intervals
#' @export
#' @importFrom stats ARMAacf
excess_stats <- function(counts,
                         intervals,
                         expected = NULL,
                         exclude = NULL,
                         model = c("quasipoisson", "poisson", "correlated"),
                         trend.nknots = 1/5,
                         harmonics = 2,
                         day.effect = TRUE,
                         control.dates = NULL,
                         max.control = 5000,
                         order.max = 30,
                         aic = TRUE,
                         verbose = TRUE){

  correlated.errors <- match.arg(model) == "correlated"

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

  ## if expected not supplied compute it
  if(is.null(expected)){
    expected <- compute_expected(counts, exclude = exclude,
                                 trend.nknots = trend.nknots,
                                 harmonics = harmonics,
                                 day.effect = day.effect)
  }
  ## compute_expected always uses quasipoisson
  if(match.arg(model) == "poisson") expected$dispersion <- 1

  if(correlated.errors){
    arfit <- fit_ar(expected, control.dates, order.max = order.max, aic = aic)
    s <- arfit$sigma
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
      s2 <- mu*expected$dispersion
    }

    total <- sum(obs)
    expected <- sum(mu)
    if(correlated.errors){
      se <- sqrt(matrix(mu, nrow = 1) %*%  Sigma %*%  matrix(mu, ncol = 1))
    } else{
      se <- sqrt(sum(s2))
    }

    return(data.frame(start = date[1], end = date[length(ind)],
                      observed = total, expected = expected, se = se,
                      observed_death_rate = total / sum(pop) * TT * 1000,
                      expected_death_rate = expected / sum(pop) * TT * 1000,
                      se_death_rate = se / sum(pop) * TT * 1000))
  })

  res <- do.call(rbind, res)

  return(res)
}
