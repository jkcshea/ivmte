context("Test of case involving only covariates, no splines.")

library(data.table)
set.seed(10L)

##------------------------
## Run MST estimator
##------------------------

ivlike <- c(ey ~ d,
            ey ~ d + x1,
            ey ~ d + x1 + x2,
            ey ~ d + x1 + x2 | x1 + x2 + z1 + z2,
            ey ~ d | factor(z2))
components <- lists.mst(d, d, c(d, x1, x2), d, d)
subsets    <- lists.mst(, , z2 %in% c(2, 3), , z2 %in% c(2, 3))

result <- mst(ivlike = ivlike,
              data = dtcf,
              components = components,
              subset = subsets,
              propensity = p,
              m0 = ~ x1 + I(x2 * u) + I(x2 * u^2),
              m1 = ~ x1 + I(x1 * x2) + u + I(x1 * u) + I(x2 * u^2),
              uname = u,
              target = "genlate",
              genlate.lb = 0.2,
              genlate.ub = 0.7,
              obseq.tol = 1.01,
              grid.Nu = 3,
              grid.Nx = 2,
              audit.Nx = 1,
              audit.Nu = 5,
              m0.inc = TRUE,
              m1.inc = TRUE,
              mte.dec = TRUE,
              treat = d)

##-------------------------
## Write function to implement test
##-------------------------

#' Calulating population mean
#'
#' Given a distribution, this function calculates the population mean
#' for each term in a formula.
#' @param formula formula, each term of which will have its mean
#'     calculated.
#' @param disribution data.table, characterizing the distribution of
#'     the variables entering into \code{formula}.
#' @param density string, name of the variable \code{data}
#'     characterizing the density.
#' @return vector, the means for each term in \code{formula}.
popmean <- function(formula, distribution, density = "f") {
    vmat  <- model.matrix(object = formula, data = distribution)
    fmat  <- as.matrix(distribution[, "f"])
    fmat  <- fmat / sum(fmat)
    means <- t(vmat) %*% fmat
    return(means)
}

#' Generate symmetric matrix
#' 
#' Function takes in a vector of values, and constructs a symmetric
#' matrix from it. Diagonals must be included. The length of the
#' vector must also be consistent with the number of "unique" entries
#' in the symmetric matrix. Note that entries are filled in along the
#' columns (i.e. equivalent to byrow = FALSE).
#' @param values vector, the values that enter into the symmetric
#'     matrix. Dimensions will be determined automatically.
#' @return matrix.
symat <- function(values) {
    k <- length(values)
    n <- 0.5 * (-1 + sqrt(1 + 8 * k))
    m <- matrix(NA, n, n)
    m[lower.tri(m, diag=TRUE)] <- values
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
   return(m)
}

#' Function to generate integral of m0 and m1
#'
#' Function carries out integral for a polynomial of degree 3.
#' @param ub scalar, upper bound of the integral.
#' @param lb scalar, lower bound of the integral.
#' @param coef vector, polynomial coefficients.
#' @return scalar.
m.int <- function(ub, lb, coef){
  return(coef[1] * (ub - lb) +
           (coef[2]/2) * (ub^2 - lb^2) +
           (coef[3]/3) * (ub^3 - lb^3))
}

#' Function to generate gamma moments
#' 
#' This function generates the gamma moments from a population level
#' data set.
#' @param data data.table.
#' @param s0 variable name (contained in the data) for the S-weight
#'     used to generate the Gamma moments for the control group.
#' @param s0 variable name (contained in the data) for the S-weight
#'     used to generate the Gamma moments for the treated group.
#' @param lb scalar, lower bound for integration.
#' @param ub scalar, upper bound for integration.
#' @return list, contains the vectors of the Gamma moments for control
#'     and treated observations.
gengamma <- function(data, s0, s1, lb, ub) {

    data <- copy(data)

    ## Gammas for D = 0
    if (!hasArg(lb) | !hasArg(ub)) {
        data[, g.0.d.0 := (1 - p) * eval(s0)]
        data[, g.0.d.1 := (1 - p) * x1 * eval(s0)]
        data[, g.0.d.2 := m.int(1, p, c(0, 1, 0)) * x2 * eval(s0)]
        data[, g.0.d.3 := m.int(1, p, c(0, 0, 1)) * x2 * eval(s0)]
    } else {
        data[, g.0.d.0 := (ub - lb) * eval(s0)]
        data[, g.0.d.1 := (ub - lb) * x1 * eval(s0)]
        data[, g.0.d.2 := m.int(ub, lb, c(0, 1, 0)) * x2 * eval(s0)]
        data[, g.0.d.3 := m.int(ub, lb, c(0, 0, 1)) * x2 * eval(s0)]
    }

    g.0.d.0 <- popmean(~ 0 + g.0.d.0, data)
    g.0.d.1 <- popmean(~ 0 + g.0.d.1, data)
    g.0.d.2 <- popmean(~ 0 + g.0.d.2, data)
    g.0.d.3 <- popmean(~ 0 + g.0.d.3, data)

    ## Gammas for D = 1
    if (!hasArg(lb) | !hasArg(ub)) {
        data[, g.1.d.0 := p * eval(s1)]
        data[, g.1.d.1 := p * x1 * eval(s1)]
        data[, g.1.d.2 := p * x1 * x2 * eval(s1)]
        data[, g.1.d.3 := m.int(p, 0, c(0, 1, 0)) * eval(s1)]
        data[, g.1.d.4 := m.int(p, 0, c(0, 1, 0)) * x1 * eval(s1)]
        data[, g.1.d.5 := m.int(p, 0, c(0, 0, 1)) * x2 * eval(s1)]
    } else {
        data[, g.1.d.0 := (ub - lb) * eval(s1)]
        data[, g.1.d.1 := (ub - lb) * x1 * eval(s1)]
        data[, g.1.d.2 := (ub - lb) * x1 * x2 * eval(s1)]
        data[, g.1.d.3 := m.int(ub, lb, c(0, 1, 0)) * eval(s1)]
        data[, g.1.d.4 := m.int(ub, lb, c(0, 1, 0)) * x1 * eval(s1)]
        data[, g.1.d.5 := m.int(ub, lb, c(0, 0, 1)) * x2 * eval(s1)]
    }

    g.1.d.0 <- popmean(~ 0 + g.1.d.0, data)
    g.1.d.1 <- popmean(~ 0 + g.1.d.1, data)
    g.1.d.2 <- popmean(~ 0 + g.1.d.2, data)
    g.1.d.3 <- popmean(~ 0 + g.1.d.3, data)
    g.1.d.4 <- popmean(~ 0 + g.1.d.4, data)
    g.1.d.5 <- popmean(~ 0 + g.1.d.5, data)

    g.0 <- c(g.0.d.0, g.0.d.1, g.0.d.2, g.0.d.3)
    g.1 <- c(g.1.d.0, g.1.d.1, g.1.d.2, g.1.d.3, g.1.d.4, g.1.d.5)

    return(list(g0 = g.0, g1 = g.1))
}

#' IV-like weighting function, OLS specification 1
#'
#' IV-like weighting function for OLS specification 1.
#' @param d 0 or 1, indicating treatment or control.
#' @return scalar.
s.ols1.d <- function(d) {
    if (d == 1) return(as.numeric(t(c(0, 1)) %*% solve(ols1.exx[1:2, 1:2])
                                  %*% c(1, 1)))
    if (d == 0) return(as.numeric(t(c(0, 1)) %*% solve(ols1.exx[1:2, 1:2])
                                  %*% c(1, 0)))
}

#' IV-like weighting function, OLS specification 2
#'
#' IV-like weighting function for OLS specification 2.
#' @param x vector, the value of the covariates other than the
#'     intercept and the treatment indicator.
#' @param d 0 or 1, indicating treatment or control.
#' @return scalar.
s.ols2.d <- function(x, d) {
    if (d == 1) return(as.numeric(t(c(0, 1, 0)) %*% solve(ols1.exx[1:3, 1:3])
                                  %*% c(1, 1, x)))
    if (d == 0) return(as.numeric(t(c(0, 1, 0)) %*% solve(ols1.exx[1:3, 1:3])
                                  %*% c(1, 0, x)))
}

#' IV-like weighting function, OLS specification 3
#'
#' IV-like weighting function for OLS specification 3.
#' @param x vector, the value of the covariates other than the
#'     intercept and the treatment indicator.
#' @param d 0 or 1, indicating treatment or control.
#' @param j scalar, position of the component one is interested in
#'     constructing the IV-like weight for.
#' @return scalar.
s.ols3 <- function(x, d, j) {
    cvec    <- replicate(4, 0)
    cvec[j] <- 1

    if (d == 1) return(as.numeric(t(cvec) %*% solve(ols2.exx) %*% c(1, 1, x)))
    if (d == 0) return(as.numeric(t(cvec) %*% solve(ols2.exx) %*% c(1, 0, x)))
}

#' IV-like weighting function, TSLS specification
#'
#' IV-like weighting function for TSLS specification.
#' @param z vector, the value of the instrument.
#' @param j scalar, position of the component one is interested in
#'     constructing the IV-like weight for.
#' @return scalar.
s.tsls <- function(z, j) {
    cvec    <- replicate(4, 0)
    cvec[j] <- 1

    return(as.numeric(t(cvec) %*%
                      solve(tsls.pi %*% t(tsls.exz)) %*%
                      tsls.pi %*% c(1, z)))
}

#' IV-like weighting function, Wald specification
#'
#' IV-like weighting function for OLS specification 2.
#' @param z vector, the value of the instrument.
#' @return scalar.
s.wald <- function(z) {
    return(((z == 3) / p.z2.3 - (z == 2) / p.z2.2) /  (ed.z2.3 - ed.z2.2))
}

##------------------------
## Implement test
##------------------------

## Construct additional variables for which we need means of
dtc$ey <- dtc$ey1 * dtc$p + dtc$ey0 * (1 - dtc$p)
dtc$eyd <- dtc$ey1 * dtc$p

varlist <- ~  eyd + ey + ey0 + ey1 + p + x1 + x2 + z1 + z2 +
    I(ey * p) + I(ey * x1) + I(ey * x2) + I(ey * z1) + I(ey * z2) +
    I(ey0 * p) + I(ey0 * x1) + I(ey0 * x2) + I(ey0 * z1) + I(ey0 * z2) +
    I(ey1 * p) + I(ey1 * x1) + I(ey1 * x2) + I(ey1 * z1) + I(ey1 * z2) +
    I(p * p) + I(p * x1) + I(p * x2) + I(p * z1) + I(p * z2) +
    I(x1 * p) + I(x1 * x1) + I(x1 * x2) + I(x1 * z1) + I(x1 * z2) +
    I(x2 * p) + I(x2 * x1) + I(x2 * x2) + I(x2 * z1) + I(x2 * z2) +
    I(z1 * p) + I(z1 * x1) + I(z1 * x2) + I(z1 * z1) + I(z1 * z2) +
    I(z2 * p) + I(z2 * x1) + I(z2 * x2) + I(z2 * z1) + I(z2 * z2)

mv1 <- popmean(varlist, dtc)
m1 <- as.list(mv1)
names(m1) <- rownames(mv1)

##-------------------------
## Construct OLS estimate 1, no controls
##-------------------------

ols1.exx <- symat(c(1, m1[["p"]], m1[["x1"]], m1[["x2"]],
                   m1[["p"]], m1[["I(p * x1)"]], m1[["I(p * x2)"]],
                   m1[["I(x1 * x1)"]], m1[["I(x1 * x2)"]],
                   m1[["I(x2 * x2)"]]))

ols1.exy <- matrix(c(m1[["ey"]], m1[["eyd"]], m1[["I(ey * x1)"]],
                     m1[["I(ey * x2)"]]))

ols1 <- (solve(ols1.exx[1:2, 1:2]) %*% ols1.exy[1:2])[2]

##-------------------------
## Construct OLS estimate 2, with controls
##-------------------------

ols2 <- (solve(ols1.exx[1:3, 1:3]) %*% ols1.exy[1:3])[2]

##-------------------------
## Construct OLS estimate 3, with controls and subsetting
##-------------------------

## mv2 <- popmean(varlist, dtc[z2 %in% c(2, 3), ])
mv2 <- popmean(varlist, subset(dtc, dtc$z2 %in% c(2, 3)))
m2 <- as.list(mv2)
names(m2) <- rownames(mv2)

ols2.exx <- symat(c(1, m2[["p"]], m2[["x1"]], m2[["x2"]],
                   m2[["p"]], m2[["I(p * x1)"]], m2[["I(p * x2)"]],
                   m2[["I(x1 * x1)"]], m2[["I(x1 * x2)"]],
                   m2[["I(x2 * x2)"]]))

ols2.exy <- matrix(c(m2[["ey"]], m2[["eyd"]],
                     m2[["I(ey * x1)"]],
                     m2[["I(ey * x2)"]]))

ols3 <- (solve(ols2.exx) %*% ols2.exy)[2:4]

##-------------------------
## Construct TSLS estimate
##-------------------------

tsls.exz <- matrix(c(1, m1[["p"]], m1[["x1"]], m1[["x2"]],
                     m1[["z1"]], m1[["I(z1 * p)"]], m1[["I(z1 * x1)"]],
                         m1[["I(z1 * x2)"]],
                     m1[["z2"]], m1[["I(z2 * p)"]], m1[["I(z2 * x1)"]],
                         m1[["I(z2 * x2)"]],
                     m1[["x1"]], m1[["I(x1 * p)"]], m1[["I(x1 * x1)"]],
                         m1[["I(x1 * x2)"]],
                     m1[["x2"]], m1[["I(x2 * p)"]], m1[["I(x2 * x1)"]],
                         m1[["I(x2 * x2)"]]),
                 nrow = 4)

tsls.ezz <- symat(c(1, m1[["z1"]], m1[["z2"]], m1[["x1"]], m1[["x2"]],
                    m1[["I(z1 * z1)"]], m1[["I(z1 * z2)"]], m1[["I(z1 * x1)"]],
                        m1[["I(z1 * x2)"]],
                    m1[["I(z2 * z2)"]], m1[["I(z2 * x1)"]], m1[["I(z2 * x2)"]],
                    m1[["I(x1 * x1)"]], m1[["I(x1 * x2)"]],
                    m1[["I(x2 * x2)"]]))

tsls.pi <- tsls.exz %*% solve(tsls.ezz)

tsls.ezy <- matrix(c(m1[["ey"]], m1[["I(ey * z1)"]], m1[["I(ey * z2)"]],
                     m1[["I(ey * x1)"]], m1[["I(ey * x2)"]]))

tsls <- (solve(tsls.pi %*% t(tsls.exz)) %*% tsls.pi %*% tsls.ezy)[2]

##-------------------------
## Construct simple Wald estimate (i.e. Wald pertaining to single instrument)
##-------------------------

## ey.z2.2 <- popmean(~ 0 + ey, dtc[z2 == 2, ])
## ey.z2.3 <- popmean(~ 0 + ey, dtc[z2 == 3, ])
ey.z2.2 <- popmean(~ 0 + ey, subset(dtc, dtc$z2 == 2))
ey.z2.3 <- popmean(~ 0 + ey, subset(dtc, dtc$z2 == 3))

## ed.z2.2 <- popmean(~ 0 + p, dtc[z2 == 2, ])
## ed.z2.3 <- popmean(~ 0 + p, dtc[z2 == 3, ])
ed.z2.2 <- popmean(~ 0 + p, subset(dtc, dtc$z2 == 2))
ed.z2.3 <- popmean(~ 0 + p, subset(dtc, dtc$z2 == 3))

wald <- (ey.z2.3 - ey.z2.2) / (ed.z2.3 - ed.z2.2)

##-------------------------
## Test equivalence of IV-like estimates
##-------------------------

test_that("IV-like estimates", {
    expect_equal(as.numeric(result$sset$s1$beta), as.numeric(ols1))
    expect_equal(as.numeric(result$sset$s2$beta), as.numeric(ols2))
    expect_equal(as.numeric(result$sset$s3$beta), as.numeric(ols3[1]))
    expect_equal(as.numeric(result$sset$s4$beta), as.numeric(ols3[2]))
    expect_equal(as.numeric(result$sset$s5$beta), as.numeric(ols3[3]))
    expect_equal(as.numeric(result$sset$s6$beta), as.numeric(tsls))
    expect_equal(as.numeric(result$sset$s7$beta), as.numeric(wald))
})

##-------------------------
## Construct gamma terms for OLS, no controls
##-------------------------

## Generate weights
dtc$s.ols1.0.d <- s.ols1.d(0)
dtc$s.ols1.1.d <- s.ols1.d(1)

## Gammas for D = 0
## m0 = ~ x1 + I(x2 * u) + I(x2 * u^2),
## ols1.0.d.0 means "OLS specification 2. For D = 0. For variable
## "d". For term 0 in md."
g.ols1 <- gengamma(dtc, quote(s.ols1.0.d), quote(s.ols1.1.d))

##-------------------------
## Construct gamma terms for OLS, with single control
##-------------------------

## Generate weights
dtc$s.ols2.0.d <- sapply(dtc$x1, s.ols2.d, d = 0)
dtc$s.ols2.1.d <- sapply(dtc$x1, s.ols2.d, d = 1)

g.ols2 <- gengamma(dtc, quote(s.ols2.0.d), quote(s.ols2.1.d))

##-------------------------
## Construct gamma terms for OLS, with controls and subset
##-------------------------

## Generate weights
dtc.x <- split(as.matrix(dtc[, .(x1, x2)]), seq(1, nrow(dtc)))

dtc$s.ols3.0.d <- unlist(lapply(dtc.x, s.ols3, d = 0, j = 2))
dtc$s.ols3.1.d <- unlist(lapply(dtc.x, s.ols3, d = 1, j = 2))

dtc$s.ols3.0.x1 <- unlist(lapply(dtc.x, s.ols3, d = 0, j = 3))
dtc$s.ols3.1.x1 <- unlist(lapply(dtc.x, s.ols3, d = 1, j = 3))

dtc$s.ols3.0.x2 <- unlist(lapply(dtc.x, s.ols3, d = 0, j = 4))
dtc$s.ols3.1.x2 <- unlist(lapply(dtc.x, s.ols3, d = 1, j = 4))

g.ols3.d  <- gengamma(dtc[z2 %in% c(2, 3), ],
                      quote(s.ols3.0.d),
                      quote(s.ols3.1.d))
g.ols3.x1 <- gengamma(dtc[z2 %in% c(2, 3), ],
                      quote(s.ols3.0.x1),
                      quote(s.ols3.1.x1))
g.ols3.x2 <- gengamma(dtc[z2 %in% c(2, 3), ],
                      quote(s.ols3.0.x2),
                      quote(s.ols3.1.x2))

## NOTE: the term for the constants are approximately 0, 1e-16.

##-------------------------
## Construct gamma terms for TSLS
##-------------------------

## Generate weights
dtc.z <- split(as.matrix(dtc[, .(z1, z2, x1, x2)]), seq(1, nrow(dtc)))
dtc$s.tsls.0.d <- unlist(lapply(dtc.z, s.tsls, j = 2))
dtc$s.tsls.1.d <- unlist(lapply(dtc.z, s.tsls, j = 2))

g.tsls <- gengamma(dtc, quote(s.tsls.0.d), quote(s.tsls.1.d))

##-------------------------
## Construct gamma terms for simple Wald (i.e. one instrument)
##-------------------------

p.z2.2 <- sum(dtc[z2 == 2, f])
p.z2.3 <- sum(dtc[z2 == 3, f])

## Generate weights
dtc$s.wald.0.d <- sapply(dtc$z2, s.wald)
dtc$s.wald.1.d <- sapply(dtc$z2, s.wald)

g.wald <- gengamma(dtc, quote(s.wald.0.d), quote(s.wald.1.d))

##-------------------------
## Construct target gammas
##-------------------------

## Generalized LATE
wald.ub <- 0.7
wald.lb <- 0.2

dtc[, w.genlate.1 := 1 / (wald.ub - wald.lb)]
dtc[, w.genlate.0 := - w.genlate.1]

g.star.genlate <- gengamma(dtc,
                           quote(w.genlate.0),
                           quote(w.genlate.1),
                           lb = wald.lb,
                           ub = wald.ub)

##-------------------------
## Test equivalence of Gamma terms
##-------------------------

test_that("Gamma moments", {
    expect_equal(as.numeric(c(result$gstar$g0, result$gstar$g1)),
                 as.numeric(unlist(g.star.genlate)))
    expect_equal(as.numeric(c(result$sset$s1$g0, result$sset$s1$g1)),
                 as.numeric(unlist(g.ols1)))
    expect_equal(as.numeric(c(result$sset$s2$g0, result$sset$s2$g1)),
                 as.numeric(unlist(g.ols2)))
    expect_equal(as.numeric(c(result$sset$s3$g0, result$sset$s3$g1)),
                 as.numeric(unlist(g.ols3.d)))
    expect_equal(as.numeric(c(result$sset$s4$g0, result$sset$s4$g1)),
                 as.numeric(unlist(g.ols3.x1)))
    expect_equal(as.numeric(c(result$sset$s5$g0, result$sset$s5$g1)),
                 as.numeric(unlist(g.ols3.x2)))
    expect_equal(as.numeric(c(result$sset$s6$g0, result$sset$s6$g1)),
                 as.numeric(unlist(g.tsls)))
    expect_equal(as.numeric(c(result$sset$s7$g0, result$sset$s7$g1)),
                 as.numeric(unlist(g.wald)))
})

##-------------------------
## Minimize observational equivalence
##-------------------------

estimates <- c(ols1, ols2, ols3, tsls, wald)

## Construct A matrix components

A <- rbind(c(g.ols1$g0, g.ols1$g1),
           c(g.ols2$g0, g.ols2$g1),
           c(g.ols3.d$g0, g.ols3.d$g1),
           c(g.ols3.x1$g0, g.ols3.x1$g1),
           c(g.ols3.x2$g0, g.ols3.x2$g1),
           c(g.tsls$g0, g.tsls$g1),
           c(g.wald$g0, g.wald$g1))

A.extra <- matrix(c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1),
                  byrow = TRUE,
                  nrow = 7)

## Construct monotonicity matrix components

grid <- matrix(c(1, 2, 0, 2), nrow = 2, byrow = TRUE)
grid <- Reduce("rbind",
               lapply(lapply(split(grid, c(1, 2)),
                             FUN = replicate,
                             n = 3,
                             simplify = FALSE),
                      Reduce, f = "rbind"))

grid <- cbind(grid, rep(c(0, 0.5, 1), times = 2))
rownames(grid) <- NULL
grid <- data.frame(grid)

colnames(grid) <- c("x1", "x2", "u")

grid <- grid[order(grid[, 3]), ] ## Match the grid ordering of the
                                 ## function, that way we can check
                                 ## the matrices defining the LP
                                 ## problem matches perfectly.

mono0 <- model.matrix(~ x1 + I(x2 * u) + I(x2 * u^2),
                      data = grid)
mono1 <- model.matrix(~ x1 + I(x1 * x2) + u + I(x1 * u) + I(x2 * u^2),
                      data = grid)

monoA0  <- mono0[c(4, 6, 3, 5), ] - mono0[c(2, 4, 1, 3), ]
monoA1  <- mono1[c(4, 6, 3, 5), ] - mono1[c(2, 4, 1, 3), ]

Azeroes <- matrix(0, ncol = 14, nrow = 4)

m0zeroes <- matrix(0, ncol = 6, nrow = 4)
m1zeroes <- matrix(0, ncol = 4, nrow = 4)

m0mono  <- cbind(Azeroes, monoA0, m0zeroes)
m1mono  <- cbind(Azeroes, m1zeroes, monoA1)
mtemono <- cbind(Azeroes, -monoA0, monoA1)

## Construct boundedness matrix components

maxy <- max(dtcf$ey)
miny <- min(dtcf$ey)

Bzeroes <- matrix(0, ncol = 14, nrow = nrow(grid))
b0zeroes <- matrix(0, ncol = ncol(mono0), nrow = nrow(grid))
b1zeroes <- matrix(0, ncol = ncol(mono1), nrow = nrow(grid))

m0bound <- cbind(Bzeroes, mono0, b1zeroes)
m1bound <- cbind(Bzeroes, b0zeroes, mono1)

## Construct full Gurobi model
model.o <- list()
model.o$obj <- c(replicate(14, 1), replicate(10, 0))
model.o$rhs <- c(estimates,
                 replicate(nrow(m0bound), miny),
                 replicate(nrow(m1bound), miny),
                 replicate(nrow(m0bound), maxy),
                 replicate(nrow(m1bound), maxy),
                 replicate(nrow(m0mono), 0),
                 replicate(nrow(m1mono), 0),
                 replicate(nrow(mtemono), 0))

model.o$sense <- c(replicate(7, "="),
                   replicate(nrow(m0bound), ">="),
                   replicate(nrow(m1bound), ">="),
                   replicate(nrow(m0bound), "<="),
                   replicate(nrow(m1bound), "<="),
                   replicate(nrow(m0mono), ">="),
                   replicate(nrow(m1mono), ">="),
                   replicate(nrow(mtemono), "<="))

model.o$A <- rbind(cbind(A.extra, A),
                   m0bound,
                   m1bound,
                   m0bound,
                   m1bound,
                   m0mono,
                   m1mono,
                   mtemono)

model.o$ub <- c(replicate(14, Inf), replicate(10, Inf))
model.o$lb <- c(replicate(14, 0), replicate(10, -Inf))

## Minimize observational equivalence deviation
model.o$modelsense <- "min"
minobseq <- gurobi::gurobi(model.o)$objbound

##-------------------------
## Obtain the bounds for generalized LATE
##-------------------------

tolerance <- 1.01
A.top <- c(replicate(14, 1), replicate(10, 0))

model.f <- list()
model.f$obj <- c(replicate(14, 0), g.star.genlate$g0, g.star.genlate$g1)
model.f$rhs <- c(tolerance * minobseq,
                 model.o$rhs)
model.f$sense <- c("<=",
                   model.o$sense)
model.f$A <- rbind(A.top,
                   model.o$A)
model.f$ub <- c(replicate(14, Inf), replicate(10, Inf))
model.f$lb <- c(replicate(14, 0), replicate(10, -Inf))

## Find bounds  with threshold
model.f$modelsense <- "min"
min_genlate <- gurobi::gurobi(model.f)
model.f$modelsense <- "max"
max_genlate <- gurobi::gurobi(model.f)

bound <- c(min_genlate$objval, max_genlate$objval)

##-------------------------
## Test equivalence of LP problem and bounds
##-------------------------

test_that("LP problem", {
    expect_equal(result$bound, bound)
    expect_equal(dim(result$lpresult$model$A),
             dim(model.f$A))
    expect_equal(as.numeric(result$lpresult$model$A),
                 as.numeric(model.f$A))
    expect_equal(as.numeric(result$lpresult$model$rhs),
                 as.numeric(model.f$rhs))
    expect_equal(as.numeric(result$lpresult$model$rhs),
                 as.numeric(model.f$rhs))
    expect_equal(result$lpresult$model$sense,
                 model.f$sense)
})
