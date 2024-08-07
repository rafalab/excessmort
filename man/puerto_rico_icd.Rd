% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/puerto_rico_icd.R
\docType{data}
\name{puerto_rico_icd}
\alias{puerto_rico_icd}
\alias{icd}
\title{Puerto Rico daily mortality by cause of death}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 262980 rows and 4 columns.
}
\usage{
data("puerto_rico_icd")
}
\description{
A data frame with Puerto Rico daily mortality counts, stratified by cause of death
from 1999 to 2020.
which includes the day hurricanes Maria made
landfall on 2017-09-20.
The outcome column includes the number of deaths for that day for that ICD-10 code.
The object `icd` is included to show the description of each ICD-10 code.
}
\keyword{datasets}
