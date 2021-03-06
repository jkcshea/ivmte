% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lp.R
\name{runGurobi}
\alias{runGurobi}
\title{Running Gurobi LP solver}
\usage{
runGurobi(lpobj, lpsolver.options)
}
\arguments{
\item{lpobj}{list of matrices and vectors defining the linear
programming problem.}

\item{lpsolver.options}{list, each item of the list should
correspond to an option specific to the LP solver selected.}
}
\value{
a list of the output from Gurobi. This includes the
    objective value, the solution vector, and the optimization status
    (status of \code{1} indicates successful optimization) .
}
\description{
This function solves the LP problem using the Gurobi package. The
object generated by \code{\link{lpSetup}} is compatible with the
\code{gurobi} function.
}
