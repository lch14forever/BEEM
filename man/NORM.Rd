% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emFunctions.r
\name{NORM}
\alias{NORM}
\title{NORM}
\usage{
NORM(
  tss,
  gradients,
  perturbInd,
  metadata,
  rmSp,
  params,
  ncpu = 10,
  norder = 3,
  scale = NA,
  smooth = FALSE,
  forceBreak = NULL
)
}
\arguments{
\item{tss}{matrix of relative abundances (proportions), variables in rows and sample in columns}

\item{gradients}{matrix of gradients of addative log (alr) transformed abundances}

\item{perturbInd}{matrix of perturbation indicator generated by preProcess}

\item{rmSp}{Speices removed for alr tranformation}

\item{params}{list(alpha=, beta=) estimated with BLASSO}

\item{ncpu}{number of CPUs used (default: 10)}

\item{norder}{order of spline basis (default: 3)}

\item{scale}{scale the median of all samples}

\item{smooth}{smooth the biomass after normalization}

\item{forceBreak}{force to break the trajectory to handle pulsed perturbation (or species invasion) (default: NULL)}

\item{meta}{metadata}
}
\description{
Function to calculate biomass for each sample
}
