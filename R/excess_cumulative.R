#' Compute cumulative excess deaths
#' @export
excess_cumulative <- function(fit, start, end){
  ind <- which(fit$date %in% seq(start, end, by = "day"))
  n <- length(ind)
  A <- matrix(1, n, n)
  A[upper.tri(A)] <- 0
  A <- sweep(A, 2, fit$expected[ind], FUN = "*")

  fhat <- matrix(fit$fitted[ind], ncol = 1)

  fit_excess <- A %*% fhat
  obs_excess <- cumsum(fit$observed[ind] - fit$expected[ind])
  fitted_se <- sqrt(diag(A %*% fit$x[ind,] %*% fit$betacov %*% t(A %*% fit$x[ind,])))
  if(fit$model == "correlated"){
    se <- sqrt(diag(A %*% fit$cov[ind, ind] %*% t(A)))
  } else{
    se <- sqrt(cumsum(fit$expected[ind]*fit$plugins$dispersion))
  }
  data.frame(date = fit$date[ind], observed = obs_excess, se = se,
             fitted = fit_excess, fitted_se = fitted_se)
}
