% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast.R
\name{evaluate.Forecast}
\alias{evaluate.Forecast}
\title{Evaluate the Forecast}
\usage{
\method{evaluate}{Forecast}(x, truth, title, color = NULL)
}
\arguments{
\item{x}{A forecast object.}

\item{truth}{A data frame containing time, true ILI.}
}
\value{
A ggplot object
}
\description{
Evaluate the forecast by  absolute error and BCI
width. Return a plot of AE and width.
}
