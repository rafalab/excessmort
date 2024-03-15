#' Diagnostic Plots for Model Fit
#' 
#' Check mean model fit via diagnostic figures of the model components
#' 
#' @param expected The output from `compute_expected` with `keep.components = TRUE`
#' @param start First day to show
#' @param end Last day to show
#' @param color Color for the expected curve
#' @param alpha alpha blending for points
#' 
#' @return A list with six ggplot objects: 
#' `population` is a time series plot of the population.
#' `seasonal` is a plot showing the estimated seasonal effect.
#' `trend` is a time series plot showing the estimated trend.
#' `weekday` is a plot of the estimated weekday effects if they were estimated.
#' `expected`is a time series plot of the expected values.
#' `residual` is a time series plot of the residuals.
#' 
#' @examples
#' 
#' library(dplyr)
#' library(lubridate)
#' library(ggplot2)
#' 
#' flu_season <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
#'
#' exclude_dates <- c(flu_season, seq(make_date(2020, 1, 1), today(), by = "day"))
#' 
#' res  <- cdc_state_counts %>%
#'  filter(state == "Massachusetts") %>%
#'  compute_expected(exclude = exclude_dates,
#'                   keep.components = TRUE)
#'                   
#' p <- expected_diagnostic(expected = res, alpha = 0.50)
#' 
#' p$population
#' p$seasonal
#' p$trend
#' p$expected
#' p$residual

#' @export
#' @import dplyr
#' @importFrom stats qnorm sd
#' @importFrom ggplot2 ggplot geom_ribbon ylab xlab geom_point coord_cartesian geom_line geom_errorbar labs
#' @importFrom lubridate yday year

expected_diagnostic <- function(expected, 
                                start  = NULL,
                                end    = NULL,
                                color  = "#D22B2B",
                                alpha  = 0.50){
  
  # -- Check that expected is a `compute_expected` object
  if(!"compute_expected" %in% class(expected) & !"excess_model" %in% class(expected)){
    stop("The expected argument needs to be the output of the computed_expected or excess_model functions.")  
  }
  
  if("excess_model" %in% class(expected)){
    if (!attr(expected, "keep.counts")) {
      stop("The object does not include counts. You need to run excess_model with keep.counts = TRUE.")
    } else{
      expected <- expected$counts
    }
  }
  
  
  # -- Check if keep.components == TRUE
  if (!attr(expected, "keep.components")) {
    if ("compute_expected" %in% class(expected)) {
      stop("The components were not found. Run the computed_expected function with keep.components = TRUE")
    } else{
      stop("The components were not found. Run the excess_model function with keep.counts = TRUE and keep.components = TRUE")
    }
  }
  
  # -- Check if the data frequency is daily or weekly
  if (!attr(expected, "frequency") %in% c(365, 52, 12)) stop("This function assumes monthly, weekly or daily data. This dataset has ", attr(expected, "frequency"), " counts per year.")
  
  # -- We use ggplot for the figures
  requireNamespace("ggplot2")
  
  # -- If NULL, set start and end dates
  if (is.null(start)) start <- min(expected$date)
  if (is.null(end))   end   <- max(expected$date)
  
  # -- Extracting data from the `expected` argument
  dat <- with(expected,
              tibble(date, 
                     outcome,
                     expected,
                     log_expected_se,
                     population = round(population), 
                     excluded))
  
  # -- Average number of deaths in non-excluded dates
  avg <- mean(filter(dat, excluded == FALSE)$outcome)
  
  # -- Seasonal data: ensuring dimensions match 
  if (attr(expected, "frequency") %in% 365) {
    seasonal_dat <- dat %>%
      filter(excluded == FALSE) %>%
      group_by(yday(date)) %>%
      summarize(avg_outcome = mean(outcome)) %>%
      ungroup() %>%
      setNames(c("t", "avg_outcome")) %>%
      filter(t != 366)
  } else{
    if (attr(expected, "frequency") %in% 52) {
      seasonal_dat <- dat %>%
        filter(excluded == FALSE) %>%
        group_by(week(date)) %>%
        summarize(avg_outcome = mean(outcome)) %>%
        ungroup() %>%
        setNames(c("t", "avg_outcome")) %>%
        filter(t != 53)
    } else{
      if (attr(expected, "frequency") %in% 12) {
        seasonal_dat <- dat %>%
          filter(excluded == FALSE) %>%
          group_by(month(date)) %>%
          summarize(avg_outcome = mean(outcome)) %>%
          ungroup() %>%
          setNames(c("t", "avg_outcome"))
      }
    }
  }
  
  # -- Seasonal viz
  p_seasonal <- seasonal_dat %>%
    mutate(s = attr(expected, "components")$seasonal$s) %>%
    ggplot(aes(t, avg_outcome)) +
    geom_point(alpha = alpha) +
    geom_line(aes(t, (s + 1) * avg), 
              color = color, 
              size  = 0.80) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Days of the year",
         y = "Counts",
         title = "Seasonal Component")
  
  # -- Extracting the trend component
  trend <- tibble(date = expected$date, 
                  trend = attr(expected, "components")$trend) %>%
    left_join(select(dat, date, population), by = "date")
  
  # -- Getting yearly average death counts
  the_freq <- attr(expected, "frequency")
  trend_obs <- dat %>%
    mutate(year = year(date)) %>%
    group_by(year) %>%
    summarize(outcome = mean(outcome)/mean(population)*1000*the_freq) %>%
    ungroup() %>%
    mutate(date = lubridate::make_date(year, 07, 01))
  
  # -- Long term trend viz
  p_trend <- ggplot() +
    geom_point(aes(date, outcome), 
               size  = 3,
               alpha = 0.80,
               data  = trend_obs) +
    geom_point(aes(date, outcome), 
               size  = 3,
               pch   = 1,
               data  = trend_obs) +
    geom_line(aes(date, trend), 
              color = color, 
              size  = 0.80,
              data  = trend) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Date",
         y = "Mortality rate (per 1,000)",
         title = "Long-term trend")
  
  # -- If weekday.effect = TRUE, then generate the figure
  if(attr(expected, "weekday.effect")){
    
    # -- Wrangling weekday data
    weekday_dat <- dat %>%
      mutate(wday = wday(date)) %>%
      left_join(attr(expected, "components")$weekday, by = c("wday" = "weekday")) %>%
      mutate(effect = avg * (1 + effect),
             wday   = wday(date, label = TRUE))
    
    # -- Wrangling average weekday data
    weekday_avg_dat <- dat %>%
      filter(excluded == FALSE) %>%
      mutate(wday = wday(date, label = TRUE)) %>%
      group_by(wday) %>%
      summarize(se      = sd(outcome) / sqrt(n()), 
                outcome = mean(outcome)) %>%
      ungroup()
    
    # -- Weekday viz
    p_weekday <- ggplot() +
      geom_errorbar(aes(wday, 
                        ymin = outcome - 2*se,
                        ymax = outcome + 2*se),
                    width = 0.01,
                    data  = weekday_avg_dat) +
      geom_point(aes(wday, outcome), 
                 shape = 18,
                 size  = 4,
                 color = "white",
                 data  = weekday_avg_dat) +
      geom_point(aes(wday, outcome), 
                 shape = 18,
                 size  = 4,
                 alpha = 0.80,
                 data  = weekday_avg_dat) +
      geom_point(aes(wday, outcome), 
                 pch  = 5,
                 size = 3,
                 data = weekday_avg_dat) +
      geom_point(aes(wday, effect),
                 color = color,
                 size  = 3,
                 data  = unique(select(weekday_dat, wday, effect))) +
      geom_point(aes(wday, effect),
                 pch  = 1,
                 size = 3,
                 data = unique(select(weekday_dat, wday, effect))) +
      labs(x = "Day of the week",
           y = "Counts",
           title = "Weekday Component")
  } else {
    p_weekday <- NULL
  }
  
  # -- Mean viz
  p_mean <- expected %>%
    filter(date >= start & date <= end) %>%
    mutate(lwr = exp(log(expected) - 2 * log_expected_se),
           upr = exp(log(expected) + 2 * log_expected_se)) %>%
    ggplot(aes(date, outcome)) +
    geom_point(alpha = alpha) +
    geom_ribbon(aes(ymin = lwr,
                    ymax = upr),
                fill  = color, alpha = 0.5) +
    geom_line(aes(date, expected), 
              color = color,
              size  = 0.80) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Date",
         y = "Counts",
         title = "Expected Mortality Counts")
  
  # -- Residual viz
  p_residual <- expected %>%
    filter(date >= start & date <= end) %>%
    mutate(difference  = outcome - expected,
           expected_se = expected * log_expected_se,
           lwr         = -2 * expected_se,
           upr         = 2 * expected_se) %>%
    ggplot(aes(date, difference)) +
    geom_point(alpha = alpha) +
    geom_ribbon(aes(ymin = lwr,
                    ymax = upr),
                fill  = color, alpha = 0.5) +
    geom_line(aes(date, 0),
              color = color,
              size  = 0.80) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Date",
         y = "Observed - expected",
         title = "Residual Plot")
  
  # -- Population plot
  p_population <- ggplot(aes(date, population), data = filter(dat, date >= start & date <= end)) +
    geom_line(size = 0.80) +
    scale_y_continuous(labels = scales::comma) +
    labs(y     = "Population",
         x     = "Date",
         title = "Population Size")
  
  # -- Add ons to poplation plot depending on the frequency
  if(attr(expected, "frequency") %in% 365) {
    p_population <- p_population +
      geom_point(aes(date, population), 
                 size  = 2,
                 shape = 21,
                 color = "black",
                 fill  = color,
                 data  = filter(dat, month(date) == 7, day(date) == 1, date >= start, date <= end))
  } else if(attr(expected, "frequency") %in% 52) {
    p_population <- p_population +
      geom_point(aes(date, population), 
                 size  = 2,
                 shape = 21,
                 color = "black",
                 fill  = color,
                 data  = filter(dat, week(date)==26, date >= start & date <= end))
  } else if(attr(expected, "frequency") %in% 12) {
    p_population <- p_population +
      geom_point(aes(date, population), 
                 size  = 2,
                 shape = 21,
                 color = "black",
                 fill  = color,
                 data  = filter(dat, month(date)==7, date >= start & date <= end))
  }
  
  # -- To return
  return(list(population = p_population,
              seasonal   = p_seasonal,
              trend      = p_trend,
              weekday    = p_weekday,
              expected   = p_mean,
              residual   = p_residual))
}