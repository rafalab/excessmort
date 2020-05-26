## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(digits = 3)

## -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
dslabs::ds_theme_set()
library(excessmort)

## -----------------------------------------------------------------------------
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
unique(cook_demographics$agegroup)

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
counts <- cdc_state_counts %>% filter(state == "Massachusetts") %>%
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
fit <- cdc_state_counts %>% filter(state == "Massachusetts") %>%
  excess_model(exclude = exclude_dates,
               start = min(.$date),
               end = max(.$date),
               knots.per.year = 12,
               weekday.effect = FALSE,
               verbose = FALSE)

## -----------------------------------------------------------------------------
excess_plot(fit)

## -----------------------------------------------------------------------------
fit$detected_intervals

## -----------------------------------------------------------------------------
cumu <- excess_cumulative(fit, 
                          start = make_date(2020, 3, 1),
                          end = make_date(2020, 5, 2))
cumu %>%
  ggplot(aes(date)) +
  geom_ribbon(aes(ymin = observed- 2*sd, ymax = observed + 2*sd), alpha = 0.5) +
  geom_point(aes(y = observed)) 


## -----------------------------------------------------------------------------
intervals <- list(flu = seq(make_date(2017, 12, 16), make_date(2018, 2, 10), by = "day"),
                  covid19 = seq(make_date(2020, 3, 14), max(cdc_state_counts$date), by = "day"))
cdc_state_counts %>% filter(state == "Massachusetts") %>%
  excess_model(exclude = exclude_dates,
               interval = intervals,
               weekday.effect = FALSE,
               verbose = FALSE)

## -----------------------------------------------------------------------------
data("puerto_rico_counts")
head(puerto_rico_counts)

## -----------------------------------------------------------------------------
counts <- collapse_counts_by_age(puerto_rico_counts, 
                                 breaks = c(0, 5, 20, 40, 60, 75, Inf))

## -----------------------------------------------------------------------------
counts <- filter(counts, agegroup == "75-Inf")

## -----------------------------------------------------------------------------
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

## -----------------------------------------------------------------------------
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")

## -----------------------------------------------------------------------------
interval_start <- c(hurricane_dates[2],
                  hurricane_dates[3],
                  Chikungunya = make_date(2014, 8, 1),
                  Covid_19 = make_date(2020, 1, 1))
before <-c(365, 365, 365, 548) ##days before event to include
after <-c(365, 365, 365, 90) # days after event to incldue

## -----------------------------------------------------------------------------
disc <- c(TRUE, TRUE, FALSE, FALSE)

## -----------------------------------------------------------------------------
f <- lapply(seq_along(interval_start), function(i){
  excess_model(counts,
               event = interval_start[i],
               start = interval_start[i] - before[i],
               end = interval_start[i] + after[i],
               exclude = exclude_dates,
               control.dates = control_dates,
               knots.per.year = 12,
               discontinuity = disc[i],
               model = "correlated")
})

## -----------------------------------------------------------------------------
excess_plot(f[[1]], title = names(interval_start)[1])

## -----------------------------------------------------------------------------
excess_plot(f[[2]],  title = names(interval_start)[2])

## -----------------------------------------------------------------------------
excess_plot(f[[3]],  title = names(interval_start)[3])

## -----------------------------------------------------------------------------
excess_plot(f[[4]],  title = names(interval_start)[4])

## -----------------------------------------------------------------------------
ndays <- 365 ## days after event to include
cumu <- lapply(seq_along(interval_start), function(i){
      excess_cumulative(f[[i]],
                      start = interval_start[i],
                      end = pmin(make_date(2020, 3, 31), interval_start[i] + ndays)) %>%
      mutate(event_day = interval_start[i], event = names(interval_start)[i])
})
cumu <- do.call(rbind, cumu)
cumu %>%
  mutate(day = as.numeric(date - event_day)) %>%
  ggplot(aes(color = event, fill = event)) +
  geom_ribbon(aes(day, ymin = fitted - 2*se, ymax = fitted + 2*se), alpha = 0.25) +
  geom_point(aes(day, observed), alpha = 0.25, cex = 1) +
  geom_line(aes(day, fitted)) +
  ggtitle("Cumulative Excess Mortality")

