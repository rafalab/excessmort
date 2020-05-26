## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(digits = 3)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(dplyr)
library(ggplot2)
library(excessmort)
dslabs::ds_theme_set()

data("cook_records")
head(cook_records)

## -----------------------------------------------------------------------------
head(cook_demographics)

## -----------------------------------------------------------------------------
counts <- compute_counts(cook_records)
head(counts)

## -----------------------------------------------------------------------------
counts <- compute_counts(cook_records, demo = cook_demographics)
head(counts)

## -----------------------------------------------------------------------------
counts <- compute_counts(cook_records, by = "agegroup", demo = cook_demographics, 
                         breaks = c(0,20,40,60,80,Inf))
head(counts)

## -----------------------------------------------------------------------------
counts <- compute_counts(cook_records, by = c("agegroup", "race", "sex"), 
                         demo = cook_demographics, 
                         breaks = c(0,20,40,60,80,Inf))
head(counts)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(dplyr)
library(lubridate)

exclude_dates <- c(seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day"),
                   seq(make_date(2020, 1, 1), max(cdc_state_counts$date), by = "day"))

## -----------------------------------------------------------------------------
counts <- cdc_state_counts %>% 
  filter(state == "Massachusetts") %>%
  compute_expected(exclude = exclude_dates, weekday.effect = FALSE)

head(counts)

## -----------------------------------------------------------------------------
expected_plot(counts)

## -----------------------------------------------------------------------------
attr(counts, "dispersion")

## -----------------------------------------------------------------------------
res  <- cdc_state_counts %>% filter(state == "Massachusetts") %>%
  compute_expected(exclude = exclude_dates, weekday.effect = FALSE,
                   keep.components = TRUE)

## -----------------------------------------------------------------------------
qplot(res$counts$date, res$trend, geom = "line", 
      xlab = "Date", ylab = "Death rate")

## -----------------------------------------------------------------------------
qplot(day, s, data = res$seasonal, geom = "line", 
      xlab = "Day of the year", ylab = "Seconal effect")

## -----------------------------------------------------------------------------
fit <- cdc_state_counts %>% 
  filter(state == "Massachusetts") %>%
  excess_model(exclude        = exclude_dates,
               start          = min(.$date),
               end            = max(.$date),
               knots.per.year = 12,
               weekday.effect = FALSE,
               verbose        = FALSE)

## -----------------------------------------------------------------------------
excess_plot(fit)

## -----------------------------------------------------------------------------
fit$detected_intervals

## -----------------------------------------------------------------------------
cumulative_deaths  <- excess_cumulative(fit, 
                                        start = make_date(2020, 03, 01),
                                        end   = make_date(2020, 05, 02))
cumulative_deaths %>%
  ggplot(aes(date)) +
  geom_ribbon(aes(ymin = observed- 2*sd, ymax = observed + 2*sd), alpha = 0.5) +
  geom_point(aes(y = observed), size=1) +
  geom_line(aes(y = observed))

## -----------------------------------------------------------------------------
intervals <- list(flu     = seq(make_date(2017, 12, 16), make_date(2018, 2, 10), by = "day"),
                  covid19 = seq(make_date(2020, 03, 14), max(cdc_state_counts$date), by = "day"))

cdc_state_counts %>% 
  filter(state == "Massachusetts") %>%
  excess_model(exclude        = exclude_dates,
               interval       = intervals,
               weekday.effect = FALSE,
               verbose        = FALSE)

