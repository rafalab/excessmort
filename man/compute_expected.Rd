% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_expected.R
\name{compute_expected}
\alias{compute_expected}
\title{Compute expected counts for each day}
\usage{
compute_expected(
  counts,
  exclude = NULL,
  include.trend = TRUE,
  trend.knots.per.year = 1/7,
  extrapolate = TRUE,
  harmonics = 2,
  frequency = NULL,
  weekday.effect = FALSE,
  keep.components = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{counts}{A data frame with dates, counts, and population size.}

\item{exclude}{A list of dates to exclude when fitting the model. This is typically the period for which you will later estimate excess counts.}

\item{include.trend}{Logical that determines if a slow trend is included in the model.}

\item{trend.knots.per.year}{Number of knots per year used for the time trend}

\item{extrapolate}{Logical that determines if the slow trend is extrapolated past the range of data used to fit. This}

\item{harmonics}{Number of harmonics to include in the seasonal effect}

\item{frequency}{Number of data points per year. If not provided, the function attempts to estimate it}

\item{weekday.effect}{A logical that determines if a day of the week effect is included in the model}

\item{keep.components}{A logical that if `TRUE` forces the function to return the estimated trend, seasonal model, and weekday effect, if included in the model.}

\item{verbose}{A logical that if `TRUE` makes function prints out updates on the estimation procedure}
}
\value{
The `counts` data.frame with two columns added: `expected` and `excluded`. 
The `expected` column is the estimated expected value of the counts for that date.
The `excluded` column is a logical vector denoting if that date was excluded when
estimating the expected value.

If the argument `keep.components` is `TRUE` a list is returned with `counts`
data.frame in the first component, the estimated trend in the second, the 
estimated seasonal effect in the third and the estimated weekday effects in the fourth.
}
\description{
Compute the expected death count for each unit of time. 
We assume counts are over-dispersed Poisson distributed with a 
trend that accounts for slow, year-to-year, changes in death rate across time and a seasonal effect. 
The function takes a data frame with 
dates and counts and returns the data frame with the expected counts as a new 
column. It also returns a logical column that is `TRUE` if that entry was 
used in the estimation procedure.
}
\details{
Periods for which excess deaths will be estimated should be excluded when estimating expected counts. These can be supplied via the `exclude` argument. 
Note that If `extrapolate` is `TRUE`, the default, the time trend will be extrapolated following the estimated trend. 
If `extrapolate` is `FALSE` the trend is assumed to be a constant equal to the estimate on the last day before extrapolation.
If the period for which excess deaths are estimated is long, extrapolation should be used with caution. 
We highly recommend exploring the estimated expected counts with the `expected_diagnostics` function.
}
\examples{
data(new_jersey_counts)
exclude_dates <- as.Date("2012-10-29") + 0:180
counts <- compute_expected(new_jersey_counts, exclude = exclude_dates, weekday.effect = TRUE)
library(ggplot2)
expected_plot(counts)

}
