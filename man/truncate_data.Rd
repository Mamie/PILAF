% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BNPR_forecast.R
\name{truncate_data}
\alias{truncate_data}
\title{Truncates a phylogentic tree}
\usage{
truncate_data(tree, truncation_time, include_trunc = F, n = 2,
  thresh = 0.01)
}
\arguments{
\item{tree}{a phylo object}

\item{truncation_time}{The time of last sampling point}

\item{include_trunc}{If to include the truncation time in the data}

\item{n}{The number of samples near the truncation time}

\item{thresh}{Maximum distance from truncation time of the n samples}
}
\description{
Truncates a phylogentic tree
}
