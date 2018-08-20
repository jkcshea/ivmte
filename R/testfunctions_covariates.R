#' Calulating population mean
#'
#' Given a distribution, this function calculates the population mean
#' for each term in a formula.
#' @param formula formula, each term of which will have its mean
#'     calculated.
#' @param distribution data.table, characterizing the distribution of
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
#' @param s1 variable name (contained in the data) for the S-weight
#'     used to generate the Gamma moments for the treated group.
#' @param lb scalar, lower bound for integration.
#' @param ub scalar, upper bound for integration.
#' @return list, contains the vectors of the Gamma moments for control
#'     and treated observations.
genGamma <- function(data, s0, s1, lb, ub) {

    ## Gammas for D = 0
    if (!hasArg(lb) | !hasArg(ub)) {
        data$g.0.d.0 <- (1 - data$p) * data[, s0]
        data$g.0.d.1 <- (1 - data$p) * data$x1 * data[, s0]
        data$g.0.d.2 <- m.int(1, data$p, c(0, 1, 0)) * data$x2 * data[, s0]
        data$g.0.d.3 <- m.int(1, data$p, c(0, 0, 1)) * data$x2 * data[, s0]
    } else {
        data$g.0.d.0 <- (ub - lb) * data[, s0]
        data$g.0.d.1 <- (ub - lb) * data$x1 * data[, s0]
        data$g.0.d.2 <- m.int(ub, lb, c(0, 1, 0)) * data$x2 * data[, s0]
        data$g.0.d.3 <- m.int(ub, lb, c(0, 0, 1)) * data$x2 * data[, s0]        
    }

    g.0.d.0 <- popmean(~ 0 + g.0.d.0, data)
    g.0.d.1 <- popmean(~ 0 + g.0.d.1, data)
    g.0.d.2 <- popmean(~ 0 + g.0.d.2, data)
    g.0.d.3 <- popmean(~ 0 + g.0.d.3, data)

    ## Gammas for D = 1
    if (!hasArg(lb) | !hasArg(ub)) {
        data$g.1.d.0 <- data$p * data[, s1]
        data$g.1.d.1 <- data$p * data$x1 * data[, s1]
        data$g.1.d.2 <- data$p * data$x1 * data$x2 * data[, s1]
        data$g.1.d.3 <- m.int(data$p, 0, c(0, 1, 0)) * data[, s1]
        data$g.1.d.4 <- m.int(data$p, 0, c(0, 1, 0)) * data$x1 * data[, s1]
        data$g.1.d.5 <- m.int(data$p, 0, c(0, 0, 1)) * data$x2 * data[, s1]
    } else {
        data$g.1.d.0 <- (ub - lb) * data[, s1]
        data$g.1.d.1 <- (ub - lb) * data$x1 * data[, s1]
        data$g.1.d.2 <- (ub - lb) * data$x1 * data$x2 * data[, s1]
        data$g.1.d.3 <- m.int(ub, lb, c(0, 1, 0)) * data[, s1]
        data$g.1.d.4 <- m.int(ub, lb, c(0, 1, 0)) * data$x1 * data[, s1]
        data$g.1.d.5 <- m.int(ub, lb, c(0, 0, 1)) * data$x2 * data[, s1]
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
#' @param exx the matrix E[XX']
#' @return scalar.
s.ols1.d <- function(d, exx) {
    if (d == 1) return(as.numeric(t(c(0, 1)) %*% solve(exx[1:2, 1:2])
                                  %*% c(1, 1)))
    if (d == 0) return(as.numeric(t(c(0, 1)) %*% solve(exx[1:2, 1:2])
                                  %*% c(1, 0)))
}

#' IV-like weighting function, OLS specification 2
#'
#' IV-like weighting function for OLS specification 2.
#' @param x vector, the value of the covariates other than the
#'     intercept and the treatment indicator.
#' @param d 0 or 1, indicating treatment or control.
#' @param exx the matrix E[XX']
#' @return scalar.
s.ols2.d <- function(x, d, exx) {
    if (d == 1) return(as.numeric(t(c(0, 1, 0)) %*% solve(exx[1:3, 1:3])
                                  %*% c(1, 1, x)))
    if (d == 0) return(as.numeric(t(c(0, 1, 0)) %*% solve(exx[1:3, 1:3])
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
#' @param exx the matrix E[XX']
#' @return scalar.
s.ols3 <- function(x, d, j, exx) {
    cvec    <- replicate(4, 0)
    cvec[j] <- 1

    if (d == 1) return(as.numeric(t(cvec) %*% solve(exx) %*% c(1, 1, x)))
    if (d == 0) return(as.numeric(t(cvec) %*% solve(exx) %*% c(1, 0, x)))
}

#' IV-like weighting function, TSLS specification
#'
#' IV-like weighting function for TSLS specification.
#' @param z vector, the value of the instrument.
#' @param j scalar, position of the component one is interested in
#'     constructing the IV-like weight for.
#' @param exz the matrix E[XZ']
#' @param pi the matrix E[XZ']E[ZZ']^{-1}
#' @return scalar.
s.tsls <- function(z, j, exz, pi) {
    cvec    <- replicate(4, 0)
    cvec[j] <- 1

    return(as.numeric(t(cvec) %*%
                      solve(pi %*% t(exz)) %*%
                      pi %*% c(1, z)))
}

#' IV-like weighting function, Wald specification
#'
#' IV-like weighting function for OLS specification 2.
#' @param z vector, the value of the instrument.
#' @param p.to P[Z = z'], where z' is value of the instrument the
#'     agent is switching to.
#' @param p.from P[Z = z], where z is the value of the instrument the
#'     agent is switching from.
#' @param e.to E[D | Z = z'], where z' is the value of the instrument
#'     the agent is switching to.
#' @param e.from E[D | Z = z], where z is the value of the instrument the
#'     agent is switching from.
#' @return scalar.
s.wald <- function(z, p.to, p.from, e.to, e.from) {
    return(((z == 3) / p.to - (z == 2) / p.from) /  (e.to - e.from))
}
