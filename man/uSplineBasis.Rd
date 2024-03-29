% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtr.R
\name{uSplineBasis}
\alias{uSplineBasis}
\title{Spline basis function}
\usage{
uSplineBasis(x, knots, degree = 0, intercept = TRUE)
}
\arguments{
\item{x}{the points to evaluate the integral of the the splines.}

\item{knots}{the knots of the spline.}

\item{degree}{the degree of the spline; default is set to 0
(constant splines).}

\item{intercept}{boolean, set to TRUE if intercept term is to be
included (i.e. an additional basis such that the sum of the
splines at every point in \code{x} is equal to 1).}
}
\value{
a matrix, the values of the integrated splines. Each row
    corresponds to a value of \code{x}; each column corresponds to
    a basis defined by the degrees and knots.
}
\description{
This function evaluates the splines that the user specifies when
declaring the MTRs. This is to be used for auditing, namely when
checking the boundedness and monotonicity conditions.
}
\examples{
## Since the splines are declared as part of the MTR, you will need
## to have parsed out the spline command. Thus, this command will be
## called via eval(parse(text = .)). In the examples below, the
## commands are parsed from the object \code{splineslist} generated
## by \code{\link[MST]{removeSplines}}. The names of the elements in
## the list are the spline commands, and the elements themselves are
## the terms that interact with the splines.

## Declare MTR function
m0 = ~ x1 + x1 : uSpline(degree = 2,
                          knots = c(0.2, 0.4)) +
    x2 : uSpline(degree = 2,
                  knots = c(0.2, 0.4)) +
    x1 : x2 : uSpline(degree = 2,
                       knots = c(0.2, 0.4)) +
    uSpline(degree = 3,
             knots = c(0.2, 0.4),
             intercept = FALSE)

## Extract spline functions from MTR function
splineslist <- removeSplines(m0)$splineslist

## Declare points at which we wish to evaluate the spline functions
x <- seq(0, 1, 0.2)

## Evaluate the splines
eval(parse(text = gsub("uSpline\\\\(",
                       "ivmte:::uSplineBasis(x = x, ",
                        names(splineslist)[1])))

eval(parse(text = gsub("uSpline\\\\(",
                       "ivmte:::uSplineBasis(x = x, ",
                       names(splineslist)[2])))
}
