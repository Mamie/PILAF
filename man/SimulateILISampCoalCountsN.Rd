% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{SimulateILISampCoalCountsN}
\alias{SimulateILISampCoalCountsN}
\title{Simulate ILI, Sampling, and Coalescent counts N Times}
\usage{
SimulateILISampCoalCountsN(n, lim, flu.Ne, flu.sampNum, ILI.sampNum,
  seed = 1)
}
\arguments{
\item{n}{Number of simulations.}

\item{lim}{A numeric tuple of start and end times.}

\item{flu.Ne}{A function of time (for flu effective popluation size).}

\item{flu.sampNum}{A numeric scalar of expected number of flu samples in lim.}

\item{ILI.sampNum}{A numeric scalar of expected number of sampling events.}
}
\description{
Simulate ILI, sampling and coalescent event counts for n Times
}
