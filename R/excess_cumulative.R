#' Compute cumulative excess deaths
#' @export
excess_cumulative <- function(fit, start, end){
  if(!"curve_fit" %in% attr(fit, "type"))
    stop("This is not the correct excess_model fit, needs curve fit.")

  ind <- which(fit$date %in% seq(start, end, by = "day"))
  n <- length(ind)
  A <- matrix(1, n, n)
  A[upper.tri(A)] <- 0
  A <- sweep(A, 2, fit$expected[ind], FUN = "*")

  fhat <- matrix(fit$fitted[ind], ncol = 1)

  fit_excess <- A %*% fhat
  obs_excess <- cumsum(fit$observed[ind] - fit$expected[ind])
  fitted_se <- sqrt(diag(A %*% fit$x[ind,] %*% fit$betacov %*% t(A %*% fit$x[ind,])))
  sd <- sqrt(diag(A %*% fit$cov[ind, ind] %*% t(A)))
  data.frame(date = fit$date[ind],
             observed = obs_excess,
             sd = sd,
             fitted = fit_excess,
             se = fitted_se)
}
