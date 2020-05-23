#' Compute year of the day ingoring Feb 29
#' 
#' @param x date
noleap_yday <- function(x){
  requireNamespace("lubridate")
  ifelse(lubridate::leap_year(x) & lubridate::month(x)>2,
         lubridate::yday(x)-1,
         lubridate::yday(x))
}
