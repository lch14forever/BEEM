% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emFunctions.r
\name{smooth.deri.pspline}
\alias{smooth.deri.pspline}
\title{smooth.deri.pspline}
\usage{
smooth.deri.pspline(x, t, deriv = 0, t.out = t, smethod = 3, norder = 3, ...)
}
\arguments{
\item{x}{vector of data}

\item{t}{time stamp of each point}

\item{deriv}{order of derivative to calculate}

\item{t.out}{time stamp of output (default: to be the same as input)}

\item{norder}{order of spline basis (default: 3)}

\item{method}{spline fitting method used for "pspline" (default: 3)}
}
\value{
a vector of fitted value after smoothing
}
\description{
Function to smooth data and calculate derivative with pspline
}
