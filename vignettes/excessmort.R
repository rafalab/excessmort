## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_options <- options()
options(digits = 3)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(knitr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(excessmort)

## ----plot-options, include=FALSE----------------------------------------------
# -- Set up for figures
theme_set(theme_bw(base_size   = 12, 
                   base_family = "Helvetica")) 

# -- Modifying plot elements globally
theme_update(
  axis.ticks        = element_line(color = "grey92"),
  axis.ticks.length = unit(.5, "lines"),
  panel.grid.minor  = element_blank(),
  legend.title      = element_text(size = 12),
  legend.text       = element_text(color = "grey30"),
  legend.background = element_rect(color = "black", fill = "white"),
  legend.key        = element_rect(fill = "white"), #FBFCFC
  legend.direction  = "horizontal",
  legend.position   = "top",
  plot.title        = element_text(size = 18, face = "bold"),
  plot.subtitle     = element_text(size = 12, color = "grey30"),
  plot.caption      = element_text(size = 9, margin = margin(t = 15)),
  plot.background   = element_rect(fill="white", color = "white"),
  panel.background  = element_rect(fill="white", color = NA),
  strip.text        = element_text(face = "bold", color = "white"),
  strip.background  = element_rect(fill = "#252525"))

## ----message=FALSE, warning=FALSE---------------------------------------------
# -- Loading Cook County records
data("cook_records")
kable(cook_records[1:6,])

## -----------------------------------------------------------------------------
# -- Cook County demographic information
kable(cook_demographics[1:6,])

## -----------------------------------------------------------------------------
# -- Aggregating death counts
counts <- compute_counts(cook_records)
kable(counts[1:6,])

## -----------------------------------------------------------------------------
# -- Aggregating death counts and computing population size from demographic data
counts <- compute_counts(cook_records, demo = cook_demographics)
kable(counts[1:6,])

## -----------------------------------------------------------------------------
# -- Aggregating death counts and computing population size by age groups
counts <- compute_counts(cook_records, by = "agegroup", demo = cook_demographics, 
                         breaks = c(0, 20, 40, 60, 80, Inf))
kable(counts[1:6,])

## -----------------------------------------------------------------------------
# -- Aggregating death counts and computing population size by age groups, race, and sex
counts <- compute_counts(cook_records, by = c("agegroup", "race", "sex"), 
                         demo = cook_demographics, 
                         breaks = c(0, 20, 40, 60, 80, Inf))
kable(counts[1:6,])

## ----message=FALSE, warning=FALSE---------------------------------------------
# -- Dates to exclude when fitting the mean model
exclude_dates <- c(seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day"),
                   seq(make_date(2020, 1, 1), max(cdc_state_counts$date), by = "day"))

## -----------------------------------------------------------------------------
# -- Fitting mean model to data from Massachusetts
counts <- cdc_state_counts %>% 
  filter(state == "Massachusetts") %>%
  compute_expected(exclude = exclude_dates)

kable(counts[1:6,])

## ----fig.align="center", fig.width=6, fig.height=4----------------------------
# -- Visualizing weekly counts and expected counts in blue
expected_plot(counts, title = "Weekly Mortality Counts in MA")

## -----------------------------------------------------------------------------
# -- Dispersion parameter from the mean model
attr(counts, "dispersion")

## ----fig.align="center", fig.width=6, fig.height=4----------------------------
# -- Fitting mean model to data from Massachusetts and retaining mean model componentss
res  <- cdc_state_counts %>% filter(state == "Massachusetts") %>%
  compute_expected(exclude = exclude_dates,
                   keep.components = TRUE)

## ----fig.align="center", fig.width=6, fig.height=4----------------------------
# -- Creating diagnostic plots
mean_diag <- expected_diagnostic(res)

# -- Trend component
mean_diag$trend

# -- Seasonal component
mean_diag$seasonal

## -----------------------------------------------------------------------------
# -- Fitting excess model to data from Massachusetts
fit <- cdc_state_counts %>% 
  filter(state == "Massachusetts") %>%
  excess_model(exclude = exclude_dates,
               start = min(.$date),
               end = max(.$date),
               knots.per.year = 12,
               verbose = FALSE)

## ----fig.align="center", fig.width=6, fig.height=4----------------------------
# -- Visualizing deviations from expected mortality in Massachusetts
excess_plot(fit, title = "Deviations from Expected Mortality in MA")

## -----------------------------------------------------------------------------
# -- Intervals of inordinate mortality found by the excess model
fit$detected_intervals

## ----fig.align="center", fig.width=6, fig.height=4----------------------------
# -- Computing excess deaths in Massachusetts from March 1, 2020 to May 9, 2020
cumulative_deaths  <- excess_cumulative(fit, 
                                        start = make_date(2020, 03, 01),
                                        end   = make_date(2020, 05, 09))

# -- Visualizing cumulative excess deaths in MA
cumulative_deaths %>%
  ggplot(aes(date)) +
  geom_ribbon(aes(ymin = observed- 2*sd, ymax = observed + 2*sd), alpha = 0.5) +
  geom_line(aes(y = observed),
            color = "white",
            size  = 1) +
  geom_line(aes(y = observed)) +
  geom_point(aes(y = observed)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x        = "Date",
       y        = "Cumulative excess deaths",
       title    = "Cumulative Excess Deaths in MA",
       subtitle = "During the first wave of Covid-19")

## -----------------------------------------------------------------------------
# -- Intervals of interest
intervals <- list(flu = seq(make_date(2017, 12, 16), make_date(2018, 2, 10), by = "day"),
                  covid19 = seq(make_date(2020, 03, 14), max(cdc_state_counts$date), by = "day"))

# -- Getting excess death statistics from the excess models for the intervals of interest
cdc_state_counts %>% 
  filter(state == "Massachusetts") %>%
  excess_model(exclude        = exclude_dates,
               interval       = intervals,
               verbose        = FALSE)

## -----------------------------------------------------------------------------
# -- Loading data from Puerto Rico
data("puerto_rico_counts")
head(puerto_rico_counts)

