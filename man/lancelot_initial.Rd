% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lancelot.R
\name{lancelot_initial}
\alias{lancelot_initial}
\title{Initial conditions for the lancelot model}
\usage{
lancelot_initial(info, n_particles, pars)
}
\arguments{
\item{info}{The result of running the \verb{$info()} method on an
initialised \link{basic} model}

\item{n_particles}{The number of particles required. Currently
only uniform initial seeding is implemented so this has no
effect}

\item{pars}{A parameter list created by \code{\link[=basic_parameters]{basic_parameters()}}; from
this list we will use the \code{population} element.}
}
\value{
A numeric vector of initial conditions
}
\description{
Create initial conditions for the lancelot model. This matches the
interface required for mcstate
}
\examples{
p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
mod <- lancelot$new(p, 0, 10)
lancelot_initial(mod$info(), 10, p)
}
