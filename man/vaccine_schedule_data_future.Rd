% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vaccination.R
\name{vaccine_schedule_data_future}
\alias{vaccine_schedule_data_future}
\title{Vaccination schedule using data and future}
\usage{
vaccine_schedule_data_future(
  data,
  region,
  uptake,
  end_date,
  mean_days_between_doses
)
}
\arguments{
\item{data}{A data.frame with columns \code{date}, \code{age_band_min},
\code{dose1} and \code{dose2}.}

\item{region}{Region to use to get total population numbers}

\item{uptake}{A vector of length 19 with fractional uptake per
group. If a single number is given it is shared across all
groups (note that this includes under-18s)}

\item{end_date}{The final day in the future to create a schedule
for. After this date the model will assume 0 vaccine doses given
so an overestimate is probably better than an underestimate.}

\item{mean_days_between_doses}{Assumed mean days between doses one
and two}
}
\value{
A \link{vaccine_schedule} object
}
\description{
Helper function to create a vaccination schedule that covers data
from the past and projects doses into the future based on the last
week of vaccination (based on JCVI order using
\link{vaccine_schedule_future}). This function is subject to change.
}
