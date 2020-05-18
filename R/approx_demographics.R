#' Interpolate demographic data
#' @export
#' @importFrom stats approx
#' @import dplyr
approx_demographics <- function(demo, first_day, last_day, by=c("year"))
{
  ################## ----- PARAMETERS ----- ##################
  # 1. demo      : Dataframe from the function getDemographics()
  # 2. first_day : First day of interpolating sequence
  # 3. last_day  : Last day of interpolating sequence

  # -- Defining xout days
  first_day <- as.Date(first_day)
  last_day  <- as.Date(last_day)
  days <- seq(first_day, last_day, by=1)

  # -- Function used within approxDemographics to interpolate population
  do_approx <- function(tab, days){

    # -- Create date variable
    tab$date <- lubridate::make_date(tab$year, 7, 1)

    # -- Variables to be return
    return(data.frame(date       = days,
                      population = approx(tab$date, tab$population, xout=days, rule=2)$y))
  }

  # -- Interpolating population
  demo <- demo %>%
    group_by_at(by) %>%
    summarize(population = sum(population, na.rm=T)) %>%
    ungroup() %>%
    group_by_at(by[-1]) %>%
    do(do_approx(., days)) %>%
    ungroup() %>%
    as_tibble()

  return(demo)
}
