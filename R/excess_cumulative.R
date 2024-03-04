#' Compute cumulative excess deaths
#' 
#' This function takes the output of the `excess_model` function, a start date, and 
#' end date and calculates excess mortality and standard errors.
#' 
#' @param fit The output of `excess_model` 
#' @param start The start date 
#' @param end The end date
#' 
#' @return A data frame with the following columns
#' \describe{
#' \item{date}{The date}
#' \item{observed}{The observed excess mortality,which is the sum of observed minus expected until that date}
#' \item{sd}{The standard deviation for excess mortality for that date if year is typical}
#' \item{fitted}{The estimated of excess mortality based on the estimated smooth event effect curve}
#' \item{se}{The standard error for `fitted`}
#' }
#'
#' @examples
#' 
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
#' excess_cumulative(f, 
#' start = as.Date("2017-12-15"), 
#' end = as.Date("2017-12-21") )
#' 
#' @export
excess_cumulative <- function(fit, start, end){
  if (!"curve_fit" %in% attr(fit, "type"))
    stop("This is not the correct excess_model fit, needs curve fit.")

  ind             <- which(fit$date %in% seq(start, end, by = "day"))
  n               <- length(ind)
  A               <- matrix(1, n, n)
  A[upper.tri(A)] <- 0
  A               <- sweep(A, 2, fit$expected[ind], FUN = "*")

  fhat <- matrix(fit$fitted[ind], ncol = 1)

  fit_excess   <- A %*% fhat
  obs_excess   <- cumsum(fit$observed[ind] - fit$expected[ind])
  varcov       <- fit$x[ind,] %*% fit$betacov %*% t(fit$x[ind,])
  diag(varcov) <- fit$se[ind]^2
  fitted_se    <- sqrt(diag(A %*%  varcov %*% t(A)))
  sd           <- sqrt(diag(A %*% (fit$cov[ind,ind] + diag(fit$log_expected_se[ind]^2)) %*% t(A)))
  data.frame(date     = fit$date[ind],
             observed = obs_excess,
             sd       = sd,
             fitted   = fit_excess,
             se       = fitted_se)
}
