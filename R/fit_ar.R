#' Fit an ar model to residuals from expected counts
#' @importFrom stats ar acf
#' @importFrom graphics plot

fit_ar <- function(expected, control.dates = NULL,
                   order.max = 5, aic = FALSE, plot = FALSE){


  if(is.null(control.dates)){
    warning("No control region suplied, using all data.")
    control.dates <- expected$date
  }

  date <- expected$date
  ind <- which(date %in% control.dates)
  mu <- expected$expected[ind]
  y <- (expected$observed[ind] - mu)/mu

  s2 <- pmax(0, mean(y^2 - 1/mu))
  w <- 1 / sqrt(1/mu + s2)

  if(order.max > 0){
    arfit <- ar(y*w, aic = aic, order.max = order.max)
  } else{
    arfit <- list(ar = c())
  }

  if(plot) acf(y*w)
  return(list(sigma = sqrt(s2), ar = arfit$ar))
}
