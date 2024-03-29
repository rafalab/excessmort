% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expected_diagnostic.R
\name{expected_diagnostic}
\alias{expected_diagnostic}
\title{Diagnostic Plots for Model Fit}
\usage{
expected_diagnostic(
  expected,
  start = NULL,
  end = NULL,
  color = "#D22B2B",
  alpha = 0.5
)
}
\arguments{
\item{expected}{The output from `compute_expected` with `keep.components = TRUE`}

\item{start}{First day to show}

\item{end}{Last day to show}

\item{color}{Color for the expected curve}

\item{alpha}{alpha blending for points}
}
\value{
A list with six ggplot objects: 
`population` is a time series plot of the population.
`seasonal` is a plot showing the estimated seasonal effect.
`trend` is a time series plot showing the estimated trend.
`weekday` is a plot of the estimated weekday effects if they were estimated.
`expected`is a time series plot of the expected values.
`residual` is a time series plot of the residuals.
}
\description{
Check mean model fit via diagnostic figures of the model components
}
\examples{

library(dplyr)
library(lubridate)
library(ggplot2)

flu_season <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")

exclude_dates <- c(flu_season, seq(make_date(2020, 1, 1), today(), by = "day"))

res  <- cdc_state_counts \%>\%
 filter(state == "Massachusetts") \%>\%
 compute_expected(exclude = exclude_dates,
                  keep.components = TRUE)
                  
p <- expected_diagnostic(expected = res, alpha = 0.50)

p$population
p$seasonal
p$trend
p$expected
p$residual
}
