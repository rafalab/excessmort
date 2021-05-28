#' Fit an ar model to residuals from expected counts
#' 
#' Helper function to estimate autoregressive mode
#' 
#' @param counts Output from `compute_excpected`
#' @param control.dates Dates to use to estimate covariance
#' @param order.max Maximum order of autoregressive process
#' @param aic Logical that determines if AIC is used
#' @param plot logical that determines if an autocorrelation plot is generated for exploration purposes
#' 
#' @importFrom stats ar acf
#' @importFrom graphics plot

fit_ar <- function(counts, control.dates = NULL,
                   order.max = 5, aic = FALSE, plot = FALSE){


  if(is.null(control.dates)){
    warning("No control region suplied, using all data.")
    control.dates <- counts$date
  }

  date       <- counts$date
  ind        <- which(date %in% control.dates)
  mu         <- counts$expected[ind]
  log_mu_var <- counts$log_expected_se[ind]^2
  y          <- (counts$outcome[ind] - mu)/mu

  s2 <- pmax(0, mean(y^2 - 1/mu - log_mu_var))
  w  <- 1 / sqrt(1/mu + s2 + log_mu_var)

  if(order.max > 0){
    arfit <- ar(y*w, aic = aic, order.max = order.max)
  } else{
    arfit <- list(ar = c())
  }

  if(plot) acf(y*w)
  return(list(sigma = sqrt(s2), ar = arfit$ar))
}
