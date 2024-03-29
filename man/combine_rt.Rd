% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support.R
\name{combine_rt}
\alias{combine_rt}
\title{Combine Rt estimates}
\usage{
combine_rt(rt, samples, rank = TRUE, weight = "infections_inc")
}
\arguments{
\item{rt}{A list of Rt calculations from
\code{\link[=lancelot_Rt_trajectories]{lancelot_Rt_trajectories()}} (though any Rt calculation that
confirms to this will work)}

\item{samples}{A list of samples from \link[mcstate:pmcmc]{mcstate::pmcmc}}

\item{rank}{A boolean deciding whether to rank trajectories by increasing
incidence or not before combining Rt estimates}

\item{weight}{A string specifying what trajectory to use for weighting,
defaults to "infections_inc" (used for Rt calculation)}
}
\value{
A list of Rt output in the same structure as the first
element of \code{rt}. All Rt estimates will be aggregated across
regions (or whatever else you are aggregating on) based on the
parameters in \code{samples}.
}
\description{
Combine Rt across multiple runs.
}
