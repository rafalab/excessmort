#' Interpolate demographic data
#' 
#' Interpolate yearly population estimates so that a population estimate is provided for 
#' each day of the year. The function `approx` is used with `rule = 2` for extrapolation.
#' 
#' @param demo A data frame with the yearly population estimates. Must have a numeric column named year and a numeric column named population.
#' @param first_day First day to interpolate. If missing the first day of the first year in demo is used.
#' @param last_day Last day to interpolate. If missing the last day of the last year in demo is used.
#' @param by Vector of column names to group by, for example different demographic strata. If missing it will extrapolate within each strata. To collapse all strata, define as by = NULL.
#' @param extrapolation.type Type of extrapolation. Either linear, constant, or none. If none is selected the NAs are returned.
#' @param ... additional parameters sent to the function `approx`.
#' @return A data frame with dates and population estimates.
#' @export
#' @importFrom stats approx
#' @import dplyr
#' 
approx_demographics <- function(demo, first_day, last_day, by,
                                extrapolation.type = c("linear", "constant", "none"), ...){
  
  extrapolation.type <- match.arg(extrapolation.type)
  
  
  if (!"year" %in% names(demo)) stop("demo is required to have a column named year.")

  if (!is.numeric(demo$year)) stop("year column must be numeric.")
  
  if (!"population" %in% names(demo)) stop("demo is required to have a column named population.")
  
  if (!is.numeric(demo$population)) stop("population column must be numeric.")
  
  if (missing(by)) by <- setdiff(colnames(demo), "population") else{ if (!"year" %in% by) by <- c("year", by) }
  
  ## If first day missing make first day of smallest year
  if (missing(first_day)) first_day <- lubridate::make_date(min(demo$year, na.rm = TRUE), 1, 1)
  
  if (missing(last_day)) last_day <- lubridate::make_date(max(demo$year, na.rm = TRUE), 12, 31)
  
  if (first_day > last_day) stop("first_day must come before last_day.")
  
  ### This is Frank Harrell's approxExtrap from the Hmisc package
  approxExtrap <- function(x, y, xout, method = "linear", n = 50, rule = 2, f = 0, ties = "ordered", na.rm = FALSE) 
  {
    if (is.list(x)) {
      y <- x[[2]]
      x <- x[[1]]
    }
    if (na.rm) {
      d <- !is.na(x + y)
      x <- x[d]
      y <- y[d]
    }
    x <- as.numeric(x)
    y <- as.numeric(y)
    d <- !duplicated(x)
    x <- x[d]
    y <- y[d]
    d <- order(x)
    x <- x[d]
    y <- y[d]
    w <- approx(x, y, xout = xout, method = method, n = n, rule = 2, f = f, ties = ties)$y
    r <- range(x)
    d <- xout < r[1]
    if (any(is.na(d))) 
      stop("NAs not allowed in xout")
    if (any(d)) 
      w[d] <- (y[2] - y[1])/(x[2] - x[1]) * (xout[d] - x[1]) + 
      y[1]
    d <- xout > r[2]
    n <- length(y)
    if (any(d)) 
      w[d] <- (y[n] - y[n - 1])/(x[n] - x[n - 1]) * (xout[d] - 
                                                       x[n - 1]) + y[n - 1]
    list(x = xout, y = w)
  }
  
  # -- Defining xout days
  first_day <- as.Date(first_day)
  last_day  <- as.Date(last_day)
  days <- seq(first_day, last_day, by = "day")
  
  # -- Function used within approxDemographics to interpolate population
  do_approx <- function(tab, days){
    
    # -- Create date variable
    tab$date <- lubridate::make_date(tab$year, 7, 1)
    
    # -- Variables to be return
    if (extrapolation.type == "linear") {
      return(data.frame(date = days,
                        population =  approxExtrap(as.numeric(tab$date), tab$population, xout = as.numeric(days), ...)$y))
    } else{
      if (extrapolation.type == "constant") {
        return(data.frame(date = days,
                          population =  approx(as.numeric(tab$date), tab$population, xout = as.numeric(days), rule = 2, ...)$y))
      } else{
        return(data.frame(date = days,
                          population =  approx(as.numeric(tab$date), tab$population, xout = as.numeric(days), rule = 1, ...)$y))
      }
    }
  }
  
  # -- Interpolating population
  demo <- demo %>%
    group_by_at(by) %>%
    summarize(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
    group_by_at(setdiff(by, "year")) %>%
    do(do_approx(., days)) %>%
    ungroup() %>%
    as_tibble()

  return(demo)
}

# Test
# demo <- expand.grid(age = factor(seq(10,90,10)), year = seq(1970, 1990, len = 21))
# demo$population <- as.numeric(demo$age)*100 + 10000 + demo$year*100
# 
# library(ggplot2)
# type <- approx_demographics(demo)
# type |> ggplot(aes(date, population, color = as.factor(age))) + geom_line() + 
#   geom_point(aes(make_date(year,7,1), population, color = as.factor(age)), data = demo)
# 
# type <- approx_demographics(demo, by = NULL)
# type |> ggplot(aes(date, population)) + geom_line()



               