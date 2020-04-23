#' Compute year of the day ingoring Feb 29
noleap_yday <- function(x){
  requireNamespace("lubridate")
  ifelse(lubridate::leap_year(x) & month(x)>2,
         lubridate::yday(x)-1,
         lubridate::yday(x))
}
