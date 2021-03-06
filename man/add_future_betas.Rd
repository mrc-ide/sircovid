% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict-support.R
\name{add_future_betas}
\alias{add_future_betas}
\alias{future_Rt}
\title{Create future betas}
\usage{
add_future_betas(sample, rt, future, rt_type = "Rt_general")

future_Rt(value, relative_to = NULL)
}
\arguments{
\item{sample}{A \code{mcstate_pmcmc} object}

\item{rt}{The results of \link{carehomes_Rt_trajectories}
with matrix elements \code{Rt_general} and \code{date} which will be used to
calculate relative Rt values.}

\item{future}{A named list of \link{future_Rt} values. Each
element name corresponds to the date that the change takes
effect (in ISO-8601 "YYYY-MM-DD" format) and must be strictly
increasing with at least two days separating changes. Note that future
values must be of length 1 or match the number of parameters in \code{rt}}

\item{rt_type}{A string giving the entry of \code{rt} to modify,
defaults to \code{Rt_general}}

\item{value}{A value or vector of values (one for each parameter set)
to add in the future. If \code{relative_to} is
\code{NULL}, then this is an \emph{absolute} Rt value, otherwise it is a
relative value.}

\item{relative_to}{Optionally an ISO 8601 (YYYY-MM-DD) format date
string or \link{Date} object indicating the date that the value
should be taken relative to.}
}
\description{
Create parameters that include beta values correspond to different
Rt scenarios in the future.
}
\details{
This function is called \code{add_future_betas} because updates beta
values, however it does this by applying Rt changes in the future
via \code{future_Rt}.

There are two basic ways that new beta values in the future might
be added. First, as relative to some previous day in the past. For
example to use a value of beta that corresponds to 1.5x the Rt
recorded on the 28th of October, we might write\preformatted{sircovid::future_Rt(1.5, "2020-10-28")
}

Second, as an absolute Rt value. For example to use a value of
beta corresponding to an Rt of 0.9 we might write\preformatted{sircovid::future_Rt(0.9)
}

These new levels need to be associated with dates at which they
occur, and will result in a piece-wise constant pattern with 1
day's transition between each level.
}
\author{
Richard Fitzjohn
}
