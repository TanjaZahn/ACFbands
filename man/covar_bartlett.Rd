% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covar_bartlett.R
\name{covar_bartlett}
\alias{covar_bartlett}
\title{Estimate Bartlett covariance matrix of autocorrelations}
\usage{
covar_bartlett(y, H, L, diag_only = FALSE)
}
\arguments{
\item{y}{a vector of observations}

\item{H}{maximum number of lags for the ACF}

\item{L}{bandwidth used for covariance estimation.}

\item{diag_only}{logical. If \code{FALSE} (default), the entire covariance matrix is estimated. Otherwise, only the diagonal elements are estimated.}
}
\value{
Either the covariance matrix of dimension \verb{H x H}(default) or vector of length \code{H}containing the diagonal elements of the covariance matrix if \code{diag_only}has been set to \code{TRUE}.
}
\description{
Estimate covariance matrix of autocorrelations according to Bartlett's formula (Bartlett, 1946) using the algorithm proposed by Mélard and Roy (1987).
}
