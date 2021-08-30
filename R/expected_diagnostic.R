#' Diagnostic Plots for Model Fit
#' 
#' Check mean model fit via diagnostic figures of the model components
#' 
#' @param expected The output from `excess_counts` with `keep.components = TRUE`
#' @param start First day to show
#' @param end Last day to show
#' @param color Color for the expected curve
#' @param alpha alpha blending for points
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
#' expected_diagnostic(expected = res, alpha = 0.50)
#' 

#' @export
#' @import dplyr
#' @importFrom stats qnorm sd
#' @importFrom ggplot2 ggplot geom_ribbon ylab xlab geom_point coord_cartesian geom_line geom_errorbar

expected_diagnostic <- function(expected, 
                                start  = NULL,
                                end    = NULL,
                                color  = "#D22B2B",
                                alpha  = 0.50){
  
  # -- Check that expected is a `compute_expected` object
  if(!"compute_expected" %in% class(expected)) stop("The expected argument needs to be the output of the computed_expected function.")
  
  # -- Check if keep.components == TRUE
  if(!"keep.components" %in% names(attributes(expected))) stop("The components were not found. Run the computed_expected function with keep.components = TRUE")
  
  # -- Check if the data frequency is daily or weekly
  if(!attr(expected$counts, "frequency") %in% c(365, 52, 12)) stop("This function assumes monthly, weekly or daily data. This dataset has ", attr(expected$counts, "frequency"), " counts per year.")
  
  # -- We use ggplot for the figures
  requireNamespace("ggplot2")
  
  # -- If NULL, set start and end dates
  if(is.null(start)) start <- min(expected$counts$date)
  if(is.null(end))   end   <- max(expected$counts$date)
  
  # -- Extracting data from the `expected` argument
  dat <- with(expected$counts,
              tibble(date, 
                     outcome,
                     expected,
                     log_expected_se,
                     population = round(population), 
                     excluded))
  
  # -- Average number of deaths in non-excluded dates
  avg <- mean(filter(dat, excluded == FALSE)$outcome)
  
  if(attr(expected$counts, "frequency") %in% 365){
    
    # -- Seasonal data
    seasonal_dat <- dat %>%
      filter(excluded == FALSE) %>%
      group_by(yday(date)) %>%
      summarize(avg_outcome = mean(outcome)) %>%
      ungroup() %>%
      setNames(c("day", "avg_outcome")) %>%
      filter(day != 366)
    
    # -- Population plot
    p_population <- ggplot(aes(date, population), data = filter(dat, date >= start & date <= end)) +
      geom_line(size = 0.80) +
      geom_point(aes(date, population), 
                 size  = 3, 
                 color = color,
                 data  = filter(dat, month(date) == 7, day(date) == 1, date >= start, date <= end)) +
      geom_point(aes(date, population), 
                 size  = 3, 
                 pch   = 1, 
                 data  = filter(dat, month(date) == 7, day(date) == 1, date >= start, date <= end)) +
      scale_y_continuous(labels = scales::comma) +
      labs(y     = "Population",
           x     = "Date",
           title = "Estimated Population Size")
    
    # -- Seasonal viz
    p_seasonal <- seasonal_dat %>%
      left_join(expected$seasonal, by = "day") %>%
      ggplot(aes(day, avg_outcome)) +
      geom_point(alpha = alpha) +
      geom_line(aes(day, (s + 1) * avg), 
                color = color, 
                size  = 0.80) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Days of the year",
           y = "Counts",
           title = "Seasonal Component")
    
    # -- Extracting the trend component
    trend <- tibble(date = expected$counts$date, trend = expected$trend) %>%
      # filter(date >= start & date <= end) %>%
      left_join(select(dat, date, population), by = "date")
    
    # -- Getting yearly average death counts
    trend_obs <- dat %>%
      # filter(date >= start & date <= end) %>%
      mutate(year = year(date)) %>%
      group_by(year) %>%
      summarize(outcome = mean(outcome)) %>%
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
      geom_line(aes(date, trend * population / (1000 * 365)), 
                color = color, 
                size  = 0.80,
                data  = trend) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Date",
           y = "Counts",
           title = "Long-term trend")
    
    # -- If weekday.effect = TRUE, then generate the figure
    if(attr(expected$counts, "weekday.effect")){
      
      # -- Wrangling weekday data
      weekday_dat <- dat %>%
        mutate(wday = wday(date)) %>%
        left_join(expected$weekday, by = c("wday" = "weekday")) %>%
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
    }
    
    # -- Mean viz
    p_mean <- expected$counts %>%
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
    p_residual <- expected$counts %>%
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
    
  } else if(attr(expected$counts, "frequency") %in% 52) {
    
    # -- Seasonal data
    seasonal_dat <- dat %>%
      filter(excluded == FALSE) %>%
      group_by(week(date)) %>%
      summarize(avg_outcome = mean(outcome)) %>%
      ungroup() %>%
      setNames(c("week", "avg_outcome"))
    
    # -- Population plot
    p_population <- ggplot(aes(date, population), data = filter(dat, date >= start & date <= end)) +
      geom_line(size = 0.80) +
      geom_point(aes(date, population), 
                 size  = 3, 
                 color = color,
                 data  = filter(dat, week(date)==26, date >= start & date <= end)) +
      geom_point(aes(date, population), 
                 size  = 3, 
                 pch   = 1, 
                 data  = filter(dat, week(date)==26, date >= start & date <= end)) +
      scale_y_continuous(labels = scales::comma) +
      labs(y     = "Population",
           x     = "Date",
           title = "Estimated Population Size")
    
    # -- Seasonal viz
    p_seasonal <- seasonal_dat %>%
      mutate(s = expected$seasonal$s) %>%
      ggplot(aes(week, avg_outcome)) +
      geom_point(alpha = alpha) +
      geom_line(aes(week, (s + 1) * avg), 
                color = color, 
                size  = 0.80) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Weeks of the year",
           y = "Counts",
           title = "Seasonal Component")
    
    # -- Extracting the trend component
    trend <- tibble(date = expected$counts$date, trend = expected$trend) %>%
      filter(date >= start & date <= end) %>%
      left_join(select(dat, date, population), by = "date")
    
    # -- Getting yearly average death counts
    trend_obs <- dat %>%
      filter(date >= start & date <= end) %>%
      mutate(year = year(date)) %>%
      group_by(year) %>%
      summarize(outcome = mean(outcome)) %>%
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
      geom_line(aes(date, trend * population / (1000 * 52)), 
                color = color, 
                size  = 0.80,
                data  = trend) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Date",
           y = "Counts",
           title = "Long-term trend")
    
    # -- Mean viz
    p_mean <- expected$counts %>%
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
    p_residual <- expected$counts %>%
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
  } else if(attr(expected$counts, "frequency") %in% 12) {
    
    # -- Seasonal data
    seasonal_dat <- dat %>%
      filter(excluded == FALSE) %>%
      group_by(month(date)) %>%
      summarize(avg_outcome = mean(outcome)) %>%
      ungroup() %>%
      setNames(c("month", "avg_outcome"))
    
    # -- Population plot
    p_population <- ggplot(aes(date, population), data = filter(dat, date >= start & date <= end)) +
      geom_line(size = 0.80) +
      geom_point(aes(date, population), 
                 size  = 3, 
                 color = color,
                 data  = filter(dat, week(date)==26, date >= start & date <= end)) +
      geom_point(aes(date, population), 
                 size  = 3, 
                 pch   = 1, 
                 data  = filter(dat, week(date)==26, date >= start & date <= end)) +
      scale_y_continuous(labels = scales::comma) +
      labs(y     = "Population",
           x     = "Date",
           title = "Estimated Population Size")
    
    # -- Seasonal viz
    p_seasonal <- seasonal_dat %>%
      mutate(s = expected$seasonal$s) %>%
      ggplot(aes(month, avg_outcome)) +
      geom_point(alpha = alpha) +
      geom_line(aes(month, (s + 1) * avg), 
                color = color, 
                size  = 0.80) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Month of the year",
           y = "Counts",
           title = "Seasonal Component")
    
    # -- Extracting the trend component
    trend <- tibble(date = expected$counts$date, trend = expected$trend) %>%
      filter(date >= start & date <= end) %>%
      left_join(select(dat, date, population), by = "date")
    
    # -- Getting yearly average death counts
    trend_obs <- dat %>%
      filter(date >= start & date <= end) %>%
      mutate(year = year(date)) %>%
      group_by(year) %>%
      summarize(outcome = mean(outcome)) %>%
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
      geom_line(aes(date, trend * population / (1000 * 12)), 
                color = color, 
                size  = 0.80,
                data  = trend) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Date",
           y = "Counts",
           title = "Long-term trend")
    
    # -- Mean viz
    p_mean <- expected$counts %>%
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
    p_residual <- expected$counts %>%
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
  }
  
  if(attr(expected$counts, "weekday.effect")){
    return(list("population" = p_population,
                "seasonal"   = p_seasonal,
                "trend"      = p_trend,
                "weekday"    = p_weekday,
                "expected"   = p_mean,
                "residual"   = p_residual))
  } else {
    return(list("population" = p_population,
                "seasonal"   = p_seasonal,
                "trend"      = p_trend,
                "expected"   = p_mean,
                "residual"   = p_residual))
  }
}