% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mst.R
\name{gmmEstimate}
\alias{gmmEstimate}
\title{Two-step GMM estimate of TE under point identification}
\usage{
gmmEstimate(sset, gstar0, gstar1, noisy = TRUE)
}
\arguments{
\item{sset}{a list of lists constructed from the function
\link{genSSet}. Each inner list should include a coefficient
corresponding to a term in an IV specification, a matrix of the
estimates of the gamma moments conditional on (X, Z) for d = 0,
and a matrix of the estimates of the gamma moments conditional
on (X, Z) for d = 1. The column means of the last two matrices
is what is used to generate the gamma moments.}

\item{gstar0}{vector, the target gamma moments for d = 0.}

\item{gstar1}{vector, the target gamma moments for d = 1.}

\item{noisy}{boolean, default set to \code{TRUE}. If \code{TRUE},
then messages are provided throughout the estimation
procedure. Set to \code{FALSE} to suppress all messages,
e.g. when performing the bootstrap.}
}
\value{
a list containing the point estimate of the treatment
    effects, the standard errors, the 90% and 95% confidence
    intervals, the convergence code (see
    \code{\link[stats]{optim}}), the coefficients on the MTR, and
    the variance/covariance matrix of the MTR coefficient
    estimates.
}
\description{
If the user sets the argument \code{point = TRUE} in the function
\code{ivmte}, then it is assumed that the treatment effect
parameter is point identified. The observational equivalence
condition is then set up as a two-step GMM problem. Solving this
GMM problem recovers the coefficients on the MTR functions m0 and
m1. Combining these coefficients with the target gamma moments
allows us to estimate the target treatment effect.
}
\examples{
## Declare empty list to be updated (in the event multiple IV like
## specifications are provided
sSet <- list()

## Declare MTR formulas
formula1 = ~ 0 + u
formula0 = ~ 0 + u

## Construct object that separates out non-spline components of MTR
## formulas from the spline components. The MTR functions are
## obtained from this object by the function 'genSSet'.
splinesList = list(removeSplines(formula0), removeSplines(formula1))

## Construct MTR polynomials
polynomials0 <- polyparse(formula = formula0,
                 data = dtm,
                 uname = u,
                 as.function = FALSE)
polynomials1 <- polyparse(formula = formula0,
                 data = dtm,
                 uname = u,
                 as.function = FALSE)

## Generate propensity score model
propensityObj <- propensity(formula = d ~ z,
                            data = dtm,
                            link = "linear")

## Generate IV estimates
ivEstimates <- ivEstimate(formula = ey ~ d | z,
                          data = dtm,
                          components = l(intercept, d),
                          treat = d,
                          list = FALSE)

## Generate target gamma moments
targetGamma <- genTarget(treat = "d",
                         m0 = ~ 1 + u,
                         m1 = ~ 1 + u,
                         uname = u,
                         target = "atu",
                         data = dtm,
                         splinesobj = splinesList,
                         pmodobj = propensityObj,
                         pm0 = polynomials0,
                         pm1 = polynomials1,
                         point = TRUE)

## Construct S-set. which contains the coefficients and weights
## corresponding to various IV-like estimands
sSet <- genSSet(data = dtm,
                sset = sSet,
                sest = ivEstimates,
                splinesobj = splinesList,
                pmodobj = propensityObj$phat,
                pm0 = polynomials0,
                pm1 = polynomials1,
                ncomponents = 2,
                scount = 1,
                yvar = "ey",
                dvar = "d",
                means = FALSE)

## Obtain point estimates using GMM
gmmEstimate(sset = sSet$sset,
            gstar0 = targetGamma$gstar0,
            gstar1 = targetGamma$gstar1)

}