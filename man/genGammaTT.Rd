% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testfunctions_covariates.R
\name{genGammaTT}
\alias{genGammaTT}
\title{Function to generate gamma moments for 'testthat'}
\usage{
genGammaTT(data, s0, s1, lb, ub)
}
\arguments{
\item{data}{data.table.}

\item{s0}{variable name (contained in the data) for the S-weight
used to generate the Gamma moments for the control group.}

\item{s1}{variable name (contained in the data) for the S-weight
used to generate the Gamma moments for the treated group.}

\item{lb}{scalar, lower bound for integration.}

\item{ub}{scalar, upper bound for integration.}
}
\value{
list, contains the vectors of the Gamma moments for control
    and treated observations.
}
\description{
This function generates the gamma moments from a population level
data set. This is specifically constructed to carry out tests.
}
