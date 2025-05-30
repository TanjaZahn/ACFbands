% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_plot.R
\name{make_plot}
\alias{make_plot}
\title{Plot autocorrelations with inference bands.}
\usage{
make_plot(rho_hat, lb, ub, color, segment)
}
\arguments{
\item{rho_hat}{A vector of length \code{H} containing the sample autocorrelations.}

\item{lb}{A vector of length \code{H} containing the lower bound of the inference band.}

\item{ub}{A vector of length \code{H} containing the upper bound of the inference band.}

\item{color}{color of the lines.}

\item{segment}{logical. Determines whether sample autocorrelations are plotted are vertically  (\code{TRUE}) or horizontally (\code{FALSE}). The first is used for significance bands, the latter for confidence bands.}
}
\value{
a vector of length \code{H} containing the estimated autocorrelations at each lag.
}
\description{
This function is used in the three main functions to create a plot.
}
