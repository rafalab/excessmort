% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/excess_cumulative.R
\name{excess_cumulative}
\alias{excess_cumulative}
\title{Compute cumulative excess deaths}
\usage{
excess_cumulative(fit, start, end)
}
\arguments{
\item{fit}{The output of `excess_model`}

\item{start}{The start date}

\item{end}{The end date}
}
\value{
A data frame with the following columns
\describe{
\item{date}{The date}
\item{observed}{The observed excess mortality,which is the sum of observed minus expected until that date}
\item{sd}{The standard deviation for excess mortality for that date if year is typical}
\item{fitted}{The estimated of excess mortality based on the estimated smooth event effect curve}
\item{se}{The standard error for `fitted`}
}
}
\description{
This function takes the output of the `excess_model` function, a start date, and 
end date and calculates excess mortality and standard errors.
}
\examples{

data(new_jersey_counts)
exclude_dates <- as.Date("2012-10-29") + 0:180
control_dates <- seq(min(new_jersey_counts$date), min(exclude_dates) - 1, by = "day")
f <- excess_model(new_jersey_counts, 
start = as.Date("2012-09-01"), 
end = as.Date("2013-09-01"),
exclude = exclude_dates,
model = "correlated",
weekday.effect = TRUE,
control.dates = control_dates)

excess_cumulative(f, 
start = as.Date("2017-12-15"), 
end = as.Date("2017-12-21") )

}
