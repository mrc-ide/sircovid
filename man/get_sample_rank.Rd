% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support.R
\name{get_sample_rank}
\alias{get_sample_rank}
\title{Get the ranking of sample by a given variable, e.g. infections}
\usage{
get_sample_rank(sample, by = "infections")
}
\arguments{
\item{sample}{An \code{mcstate_pmcmc} object}

\item{by}{The name of the variable used for ranking. The ranking is
computed based on the value of this variable at the last time step.
Default is set to "infections", which will give the rank based on the
cumulative number of infections.}
}
\value{
A vector of integers giving the ranks
}
\description{
Get the ranking of sample
}
