% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carehomes.R
\name{carehomes_compare}
\alias{carehomes_compare}
\title{Compare observed and modelled data for the carehomes model}
\usage{
carehomes_compare(state, observed, pars)
}
\arguments{
\item{state}{State vector for the end of the current day. This is
assumed to be filtered following \code{\link[=carehomes_index]{carehomes_index()}} so contains
10 rows corresponding to ICU, general beds, admissions, deaths and
seroconversion compartments.}

\item{observed}{Observed data. At the moment please see the tests
for a full list as this changes frequently (and this function
may be removed in future).}

\item{pars}{A list of parameters, as created by
\code{\link[=carehomes_parameters]{carehomes_parameters()}}}
}
\value{
A vector of log likelihoods, the same length as the number
of particles (the number of columns in the modelled state)
}
\description{
Compare observed and modelled data from the \link{carehomes} model. This
conforms to the mcstate interface.
}
