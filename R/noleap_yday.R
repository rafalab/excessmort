#' Compute year of the day ingoring Feb 29
noleap_yday <- function(x) ifelse(leap_year(x) & month(x)>2, yday(x)-1, yday(x))
