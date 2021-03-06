% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PILAF.R
\name{forecast.PILAF}
\alias{forecast.PILAF}
\title{Forecast a PILAF Object.}
\usage{
\method{forecast}{PILAF}(x, time.forecast, week.forecast,
  method = "count", formula = "NULL", return_model = F,
  verbose = F)
}
\arguments{
\item{x}{A PILAF object.}

\item{time.forecast}{A negative numeric vector of timepoints to forecast.}

\item{week.forecast}{A numeric vector of week corresponding to timepoints to forecast.}

\item{formula}{An optional character scalar of formula (debugging).}
}
\description{
Forecast a PILAF object using one of the three methods: count only, joint
modeling of count and coalescent events, joint modeling of count, coalescent
and sampling envets.
}
