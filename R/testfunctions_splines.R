#' Integrating splines
#'
#' This function simply integrates the splines.
#' @param ub scalar, upperbound of integral.
#' @param lb scalar, lowerbound of integral.
#' @param knots vector, knots of the spline.
#' @param degree scalar, degre of spline.
#' @param intercept boolean, set to TRUE if spline basis should
#'     include a component so that the basis sums to 1.
#' @return vector, each component being the integral of a basis.
splineInt <- function(ub, lb, knots, degree, intercept = FALSE) {
    splines2::ibs(x = ub,
                  knots = knots,
                  degree = degree,
                  intercept = intercept,
                  Boundary.knots = c(0, 1)) -
        splines2::ibs(x = lb,
                      knots = knots,
                      degree = degree,
                      intercept = intercept,
                      Boundary.knots = c(0, 1))
}

#' IV-like weighting function, OLS specifications
#'
#' IV-like weighting function for OLS specifications.
#' @param x vector, the value of the covariates other than the
#'     intercept and the treatment indicator.
#' @param d 0 or 1, indicating treatment or control.
#' @param j scalar, position of the component one is interested in
#'     constructing the IV-like weight for.
#' @param exx matrix corresponding to E[XX'].
#' @return scalar.
sOlsSplines <- function(x = NULL, d, j, exx) {
    if (!is.na(x)) {
        cvec    <- replicate((length(x) + 2), 0)
        cvec[j] <- 1
        if (d == 1) return(as.numeric(t(cvec) %*% solve(exx) %*% c(1, 1, x)))
        if (d == 0) return(as.numeric(t(cvec) %*% solve(exx) %*% c(1, 0, x)))
    } else {
        cvec    <- replicate(2, 0)
        cvec[j] <- 1
        if (d == 1) return(as.numeric(t(cvec) %*% solve(exx) %*% c(1, 1)))
        if (d == 0) return(as.numeric(t(cvec) %*% solve(exx) %*% c(1, 0)))
    }
}

#' IV-like weighting function, TSLS specification
#'
#' IV-like weighting function for TSLS specification.
#' @param z vector, the value of the instrument.
#' @param d 0 or 1, indicating treatment or control (redundant in this
#'     function; included to exploit apply()).
#' @param j scalar, position of the component one is interested in
#'     constructing the IV-like weight for.
#' @param exz matrix, corresponds to E[XZ'].
#' @param pi matrix, corresponds to E[XZ']E[ZZ']^{-1}, the first stage
#'     regression.
#' @return scalar.
sTslsSplines <- function(z, d, j, exz, pi) {
    cvec    <- replicate(nrow(exz), 0)
    cvec[j] <- 1
    
    return(as.numeric(t(cvec) %*%
                      solve(pi %*% t(exz)) %*% pi %*% c(1, z)))
}

#' Target weighting function, for ATT
#' 
#' Target weighting function, for the ATT.
#' @param z vector, the value of the instrument (redundant in this
#'     function; included to exploit apply()).
#' @param d 0 or 1, indicating treatment or control (redundant in this
#'     function; included to exploit apply()).
#' @param ed scalar, unconditional probability of taking up treatment.
#' @return scalar.
wAttSplines <- function(z, d, ed) {
    1 / ed
}

#' Generating the Gamma moments for splines, for 'testthat'
#'
#' This function generates the Gamma moments for a given set of
#' weights. This funciton is written specifically for tests.
#' @param distr data.frame, the distribution of the data.
#' @param weight function, the S-function corresponding to a
#'     particular IV-like estimand.
#' @param zvars vector, string names of the covariates, other than the
#'     intercept and treatment variable.
#' @param u1s1 matrix, the spline basis for the treated group ("u1")
#'     corresponding to the first (and only) spline specification
#'     ("s1").
#' @param u0s1 matrix, the spline basis for the control group ("u0")
#'     corresponding to the first spline specification ("s1").
#' @param u0s2 matrix, the spline basis for the control group ("u0")
#'     corresponding to the second spline specification ("s2").
#' @param target boolean, set to \code{TRUE} if the gamma moment being
#'     generated corresponds to the target parameter.
#' @param ... all other arguments that enter into \code{weight},
#'     excluding the argument \code{d} for treatment indicator.
#' @return vector, the Gamma moments associated with \code{weight}.
genGammaSplinesTT <- function(distr, weight, zvars, u1s1, u0s1, u0s2,
                     target = FALSE, ...) {
    if (hasArg(zvars)) {
        zmat <- as.matrix(distr[, zvars])
    } else {
        zmat <- matrix(NA, nrow = nrow(distr))
    }
    
    s0 <- apply(X = zmat, MARGIN = 1, FUN = weight, d = 0, ...)
    s1 <- apply(X = zmat, MARGIN = 1, FUN = weight, d = 1, ...)

    ## Construct m0 moments
    mu0s1 <- sweep(u0s1, MARGIN = 1, STATS = distr[, "x"], FUN = "*")
    mu0s1 <- sweep(mu0s1, MARGIN = 1, STATS = s0, FUN = "*")
    
    mu0s2 <- sweep(u0s2, MARGIN = 1, STATS = s0, FUN = "*")

    if (target == FALSE) muu2 <- (1 / 3) * (1 - distr[, "p"] ^ 3) * s0
    if (target == TRUE)  muu2 <- (1 / 3) * (distr[, "p"] ^ 3) * s0

    g0 <- c(sum(muu2 * distr[, "f"]),
            colSums(sweep(mu0s2,
                          MARGIN = 1,
                          STATS = distr[, "f"],
                          FUN = "*")),
            colSums(sweep(mu0s1,
                          MARGIN = 1,
                          STATS = distr[, "f"],
                          FUN = "*")))
    
    names(g0) <- c("I(u^2)",
                   "u0S1.1:1", "u0S1.2:1", "u0S1.3:1",
                   "u0S2.1:x", "u0S2.2:x", "u0S2.3:x", "u0S2.4:x")
    
    ## Construct m1 moments
    m1int <- s1 * distr[, "p"]
    
    m1x <- distr[, "x"] * s1 * distr[, "p"]
    
    mu1s1 <- sweep(u1s1, MARGIN = 1, STATS = s1, FUN = "*")
  
    g1 <- c(sum(m1int * distr[, "f"]),
            sum(m1x * distr[, "f"]),
            colSums(sweep(mu1s1,
                          MARGIN = 1,
                          STATS = distr[, "f"],
                          FUN = "*")))
    
    names(g1) <- c("(Intercept)", "x",
                   "u1S1.1:1", "u1S1.2:1", "u1S1.3:1", "u1S1.4:1")
    
    return(list(g0 = g0,
                g1 = g1))
}
