% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vaccination.R
\name{vaccine_schedule_from_data}
\alias{vaccine_schedule_from_data}
\title{Create historical vaccine schedule}
\usage{
vaccine_schedule_from_data(data, region, uptake, carehomes = TRUE)
}
\arguments{
\item{data}{A data.frame with columns \code{date}, \code{age_band_min},
and numbered doses columns, e.g. if there are three doses
these should be \code{dose1}, \code{dose2} and \code{dose3}. Values of
\code{age_band_min} should be either multiples of 5 or NA - the
latter means those doses are not age-specific and will be
distributed across all ages according to priority after
all age-specific doses have already been dealt with}

\item{region}{Region to use to get total population numbers}

\item{uptake}{A matrix of 19 rows, and number of columns equal to
number of doses. The (i,j)th entry gives the fractional uptake
of dose j for group i. Should be non-increasing across rows}

\item{carehomes}{Logical parameter, whether or not we have carehomes in
the model. Default is TRUE}
}
\value{
A \link{vaccine_schedule} object
}
\description{
Create a historical vaccine schedule from data
}
