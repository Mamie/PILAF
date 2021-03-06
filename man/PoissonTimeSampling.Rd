% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{PoissonTimeSampling}
\alias{PoissonTimeSampling}
\title{Poisson Process Simulation}
\usage{
PoissonTimeSampling(traj, c, lim)
}
\arguments{
\item{traj}{A function of time.}

\item{c}{A numeric as proportionality constant on traj for intensity.}

\item{lim}{A tuple that specifies the start and end time.}
}
\value{
A numeric vector of event times.
}
\description{
Simulation of Poisson process given intensity.
}
\seealso{
PILAF::ComputePropCnst
}
