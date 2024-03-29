% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vaccination.R
\name{vaccine_schedule_scenario}
\alias{vaccine_schedule_scenario}
\title{High-level vaccine scenario creation}
\usage{
vaccine_schedule_scenario(
  schedule_past,
  doses_future,
  end_date,
  mean_days_between_doses,
  priority_population,
  lag_groups = NULL,
  lag_days = NULL,
  boosters_future = NULL,
  boosters_prepend_zero = TRUE,
  booster_proportion = rep(1L, 19)
)
}
\arguments{
\item{schedule_past}{A \link{vaccine_schedule} object corresponding to
previously carried out vaccination.}

\item{doses_future}{A named vector of vaccine doses to give in the
future. Names must be in ISO date format.}

\item{end_date}{The final day in the future to create a schedule
for. After this date the model will assume 0 vaccine doses given
so an overestimate is probably better than an underestimate.}

\item{mean_days_between_doses}{Assumed mean days between doses one
and two}

\item{priority_population}{Output from
\link{vaccine_priority_population}, giving the number of people
to vaccinate in each age (row) and priority group (column)}

\item{lag_groups}{Row indices, corresponding to age
groups in which a lag should be added to the start time of the dose
schedule returned by \link{vaccine_schedule}, if NULL then no lag is added.
Ignored if \code{lag_groups} is NULL.}

\item{lag_days}{If \code{lag_groups} is not NULL then specifies the number of
days to add the start of the dose schedule for the given groups. Ignored
if \code{lag_groups} is NULL.}

\item{boosters_future}{Optional named vector of booster doses to give in
the future. Names must be in ISO date format.}

\item{boosters_prepend_zero}{If TRUE (default) and \code{boosters_future} is
not NULL then sets booster doses to zero before the first date in
\code{boosters_future}. This is in contrast to when it is FALSE and the
previous value in \code{schedule_past} is replicated until the first date in
boosters_future. Note that this should rarely be FALSE as this will
likely lead to duplicating daily doses that are already replicated in
\code{doses_future}.}

\item{booster_proportion}{Proportion of the groups in
\code{priority_population} to boost, default is all groups; ignored if
\code{booster_daily_doses_value} is NULL.}
}
\value{
A \link{vaccine_schedule} object
}
\description{
Create a vaccination scenario
}
