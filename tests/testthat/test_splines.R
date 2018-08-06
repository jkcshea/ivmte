context("Test of case involving only covariates, no splines.")

set.seed(10L)

##------------------------
## Define functions
##------------------------

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

#' Calulating population mean
#'
#' Given a distribution, this function calculates the population mean
#' for each term in a formula.
#' @param formula formula, each term of which will have its mean
#'     calculated.
#' @param disribution data.frame, characterizing the distribution of
#'     the variables entering into \code{formula}.
#' @param density string, name of the variable \code{distribution}
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
sOls <- function(x = NULL, d, j, exx) {
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
sTsls <- function(z, d, j, exz, pi) {
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
wAtt <- function(z, d, ed) {
    1 / ed
}

#' Generating the Gamma moments
#'
#' This function generates the Gamma moments for a given set of
#' weights.
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
#' @param ... all other arguments that enter into \code{weight},
#'     excluding the argument \code{d} for treatment indicator.
#' @return vector, the Gamma moments associated with \code{weight}.
genGamma <- function(distr, weight, zvars, u1s1, u0s1, u0s2,
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

##------------------------
## Run MST estimator
##------------------------

ivlike <- c(ey ~ d,
            ey ~ d + x,
            ey ~ d + x | z + x)

components <- lists.mst(c(intercept, d), d, c(d, x))
## FIX: You should address the case where an empty list is provided.

set.seed(10L)
result <- mst(ivlike = ivlike,
              data = dtsf,
              components = components,
              propensity = p,
              m1 = ~ x + uSplines(degree = 2,
                                  knots = c(0.3, 0.6),
                                  intercept = FALSE),
              m0 = ~ 0 + x : uSplines(degree = 0,
                                      knots = c(0.2, 0.5, 0.8),
                                      intercept = TRUE) +
                  uSplines(degree = 1,
                           knots = c(0.4),
                           intercept = TRUE) +
                  I(u ^ 2),
              uname = u,
              target = "att",
              obseq.tol = 1.01,
              grid.Nu = 3,
              grid.Nx = 2,
              audit.Nx = 1,
              audit.Nu = 5,
              m1.ub = 55,
              m0.lb = 0,
              mte.inc = TRUE)

## To test your custom weight, just submit a weight that matches a
## predefined one, e.g. submit a function calculating the weight for
## the ATT.

##-------------------------
## Perform tests
##-------------------------

## Construct additional variables for which we need means of
dts$ey  <- dts$ey1 * dts$p + dts$ey0 * (1 - dts$p)
dts$eyd <- dts$ey1 * dts$p

varlist <- ~  eyd + ey + ey0 + ey1 + p + x + z +
    I(ey * p) + I(ey * x) + I(ey * z) +
    I(ey0 * p) + I(ey0 * x) + I(ey0 * z) +
    I(ey1 * p) + I(ey1 * x) + I(ey1 * z) + 
    I(p * p) + I(p * x) + I(p * z) +
    I(x * p) + I(x * x) + I(x * z) +  
    I(z * p) + I(z * x) + I(z * z)

mv <- popmean(varlist, dts)
m <- as.list(mv)
names(m) <- rownames(mv)

##-------------------------
## Construct OLS estimate 1
##-------------------------

exx <- symat(c(1, m[["p"]], m[["x"]],
                    m[["p"]], m[["I(p * x)"]],
                    m[["I(x * x)"]]))

exy <- matrix(c(m[["ey"]], m[["eyd"]], m[["I(ey * x)"]]))

ols1 <- solve(exx[1:2, 1:2]) %*% exy[1:2]

##-------------------------
## Construct OLS estimate 2
##-------------------------

ols2 <- solve(exx) %*% exy

##-------------------------
## Construct TSLS estimate
##-------------------------

exz <- matrix(c(1, m[["p"]], m[["x"]],
                m[["z"]], m[["I(z * p)"]], m[["I(z * x)"]],
                m[["x"]], m[["I(x * p)"]], m[["I(x * x)"]]), 
              nrow = 3)
              
ezz <- symat(c(1, m[["z"]], m[["x"]],
               m[["I(z * z)"]], m[["I(z * x)"]],
               m[["I(x * x)"]]))

pi <- exz %*% solve(ezz)

ezy <- matrix(c(m[["ey"]], m[["I(ey * z)"]], m[["I(ey * x)"]]))

tsls <- (solve(pi %*% t(exz)) %*% pi %*% ezy)

##-------------------------
## Test equivalence of IV-like estimates
##-------------------------

test_that("IV-like estimates", {
    expect_equal(as.numeric(result$sset$s1$beta), as.numeric(ols1[1]))
    expect_equal(as.numeric(result$sset$s2$beta), as.numeric(ols1[2]))
    expect_equal(as.numeric(result$sset$s3$beta), as.numeric(ols2[2]))
    expect_equal(as.numeric(result$sset$s4$beta), as.numeric(tsls[2]))
    expect_equal(as.numeric(result$sset$s5$beta), as.numeric(tsls[3]))
})

##-------------------------
## Generate gamma terms
##-------------------------

## Spline integrals for target etimand
t1s1 <- t(sapply(X = dts$p,
                 FUN = splineInt,
                 lb = 0,
                 degree = 2,
                 knots = c(0.3, 0.6),
                 intercept = FALSE))

t0s1 <- t(sapply(X = dts$p,
                 FUN = splineInt,
                 lb = 0,
                 degree = 0,
                 knots = c(0.2, 0.5, 0.8),
                 intercept = TRUE))

t0s2 <- t(sapply(X = dts$p,
                 FUN = splineInt,
                 lb = 0,
                 degree = 1,
                 knots = c(0.4),
                 intercept = TRUE))

## Spline integrals for IV-like estimands
u1s1 <- t(sapply(X = dts$p,
                 FUN = splineInt,
                 lb = 0,
                 degree = 2,
                 knots = c(0.3, 0.6),
                 intercept = FALSE))

u0s1 <- t(sapply(X = dts$p,
                 FUN = splineInt,
                 ub = 1,
                 degree = 0,
                 knots = c(0.2, 0.5, 0.8),
                 
                 intercept = TRUE))

u0s2 <- t(sapply(X = dts$p,
                 FUN = splineInt,
                 ub = 1,
                 degree = 1,
                 knots = c(0.4),
                 intercept = TRUE))

## Target Gamma terms
gstar <- genGamma(distr = dts,
                  u1s1 = t1s1,
                  u0s1 = t0s1,
                  u0s2 = t0s2,
                  weight = wAtt,
                  ed = m[["p"]],
                  target = TRUE)
gstar$g0 <- - gstar$g0

## S-set Gamma terms

g1 <- genGamma(distr = dts,
               u1s1 = u1s1,
               u0s1 = u0s1,
               u0s2 = u0s2,
               weight = sOls,
               j = 1,
               exx = exx[1:2, 1:2])

g2 <- genGamma(distr = dts,
               u1s1 = u1s1,
               u0s1 = u0s1,
               u0s2 = u0s2,
               weight = sOls,
               j = 2,
               exx = exx[1:2, 1:2])

g3 <- genGamma(distr = dts,
               u1s1 = u1s1,
               u0s1 = u0s1,
               u0s2 = u0s2,
               zvars = "x",
               weight = sOls,
               j = 2,
               exx = exx)

g4 <- genGamma(distr = dts,
               u1s1 = u1s1,
               u0s1 = u0s1,
               u0s2 = u0s2,
               zvars = c("z", "x"),
               weight = sTsls,
               j = 2,
               exz = exz,
               pi = pi)

g5 <- genGamma(distr = dts,
               u1s1 = u1s1,
               u0s1 = u0s1,
               u0s2 = u0s2,
               zvars = c("z", "x"),
               weight = sTsls,
               j = 3,
               exz = exz,
               pi = pi)

##-------------------------
## Test equivalence of Gamma moments
##-------------------------

test_that("Gamma moments", {
    expect_equal(result$gstar, gstar)
    expect_equal(result$sset$s1[2:3], g1)
    expect_equal(result$sset$s1[2:3], g1)
    expect_equal(result$sset$s2[2:3], g2)
    expect_equal(result$sset$s3[2:3], g3)
    expect_equal(result$sset$s4[2:3], g4)
    expect_equal(result$sset$s5[2:3], g5)
})

##-------------------------
## Minimize observational equivalence
##-------------------------

estimates <- c(ols1, ols2[2], tsls[c(2, 3)])

## Construct A matrix components

A <- rbind(c(g1$g0, g1$g1),
           c(g2$g0, g2$g1),
           c(g3$g0, g3$g1),
           c(g4$g0, g4$g1),
           c(g5$g0, g5$g1))

Aextra <- matrix(0, nrow = nrow(A), ncol = 2 * nrow(A))
for (i in 1:nrow(A)) {
    Aextra[i, (i * 2 - 1)] <- -1
    Aextra[i, (i * 2)] <- 1
}

## Construct monotonicity matrix components

uGrid <- unique(result$auditgrid$grid[, "u"])
xMat <- cbind(1, rep(c(-1, 0, 1), each = length(uGrid)))
u1s1 <- splines2::bSpline(x = uGrid,
                          degree = 2,
                          knots = c(0.3, 0.6),
                          intercept = FALSE)
u1s1 <- do.call("rbind", rep(list(u1s1), 3))

u0s1 <- splines2::bSpline(x = uGrid,
                          degree = 0,
                          knots = c(0.2, 0.5, 0.8),
                          intercept = TRUE)
u0s1 <- do.call("rbind", rep(list(u0s1), 3))
u0s1 <- sweep(x = u0s1,
              MARGIN = 1,
              STATS = xMat[, 2],
              FUN = "*")

u0s2 <- splines2::bSpline(x = uGrid,
                          degree = 1,
                          knots = c(0.4),
                          intercept = TRUE)
u0s2 <- do.call("rbind", rep(list(u0s2), 3))

u2Mat <- rep(uGrid ^ 2, 3)

mono0 <- cbind(u2Mat, u0s2, u0s1)
mono1 <- cbind(xMat, u1s1)

colnames(mono0) <- c("I(u^2)",
                     "u0S1.1:1", "u0S1.2:1", "u0S1.3:1",
                     "u0S2.1:x", "u0S2.2:x", "u0S2.3:x", "u0S2.4:x")

colnames(mono1) <- c("(Intercept)", "x",
                  "u1S1.1:1", "u1S1.2:1", "u1S1.3:1", "u1S1.4:1")

monoA0 <- mono0[seq(1, 3 * length(uGrid))[-seq(1,
                                               3 * length(uGrid),
                                               by = length(uGrid))], ] -
    mono0[seq(1, 3 * length(uGrid))[-seq(length(uGrid),
                                         3 * length(uGrid),
                                         by = length(uGrid))], ]
monoA1 <- mono1[seq(1, 3 * length(uGrid))[-seq(1,
                                               3 * length(uGrid),
                                               by = length(uGrid))], ] -
    mono1[seq(1, 3 * length(uGrid))[-seq(length(uGrid),
                                         3 * length(uGrid),
                                         by = length(uGrid))], ]

Azeroes <- matrix(0, ncol = ncol(Aextra), nrow = nrow(monoA0))
mtemono <- cbind(Azeroes, -monoA0, monoA1)

## Construct boundedness matrix components

maxy <- max(c(dts$ey0, dts$ey1))
miny <- min(c(dts$ey0, dts$ey1))

Bzeroes <- matrix(0, ncol = ncol(Aextra), 3 * length(uGrid))

b0zeroes <- matrix(0, ncol = ncol(mono0), nrow = 3 * length(uGrid))
b1zeroes <- matrix(0, ncol = ncol(mono1), nrow = 3 * length(uGrid))

## Match the order of the matrix generated by the function
bOrder <- c(sapply(seq(0, length(uGrid) - 1),
                   function(x) x + seq(2 * length(uGrid) + 1,
                                         1,                                       ,
                                         by = -length(uGrid))))

m0bound <- cbind(Bzeroes, mono0[bOrder, ], b1zeroes)
m1bound <- cbind(Bzeroes, b0zeroes, mono1[bOrder, ])

## Construct full Gurobi model
modelO <- list()
modelO$obj <- c(replicate(ncol(Aextra), 1),
                 replicate(ncol(monoA0) + ncol(monoA1), 0))

modelO$rhs <- c(estimates,
                 replicate(nrow(m0bound), 0),
                 replicate(nrow(m1bound), miny),
                 replicate(nrow(m0bound), maxy),
                 replicate(nrow(m1bound), 55),
                 replicate(nrow(mtemono), 0))

modelO$sense <- c(replicate(length(estimates), "="),
                   replicate(nrow(m0bound), ">="),
                   replicate(nrow(m1bound), ">="),
                   replicate(nrow(m0bound), "<="),
                   replicate(nrow(m1bound), "<="),
                   replicate(nrow(mtemono), ">="))

modelO$A <- rbind(cbind(Aextra, A),
                   m0bound,
                   m1bound,
                   m0bound,
                   m1bound,
                   mtemono)

modelO$ub <- c(replicate(ncol(Aextra), Inf),
                replicate(ncol(monoA0) + ncol(monoA1), Inf))
modelO$lb <- c(replicate(ncol(Aextra), 0),
                replicate(ncol(monoA0) + ncol(monoA1), -Inf))

## Minimize observational equivalence deviation
modelO$modelsense <- "min"
minobseq <- gurobi::gurobi(modelO)$objbound

##-------------------------
## Obtain the bounds for generalized LATE
##-------------------------

tolerance <- 1.01

Atop <- c(replicate(ncol(Aextra), 1),
          replicate(ncol(monoA0) + ncol(monoA1), 0))

modelF <- list()
modelF$obj <- c(replicate(ncol(Aextra), 0),
                gstar$g0,
                gstar$g1)
modelF$rhs <- c(tolerance * minobseq,
                 modelO$rhs)
modelF$sense <- c("<=",
                   modelO$sense)
modelF$A <- rbind(Atop,
                  modelO$A)
modelF$ub <- c(replicate(ncol(Aextra), Inf),
               replicate(ncol(mono0) + ncol(mono1), Inf))
modelF$lb <- c(replicate(ncol(Aextra), 0),
               replicate(ncol(mono0) + ncol(mono1), -Inf))

## Find bounds with threshold
modelF$modelsense <- "min"
min_genlate <- gurobi::gurobi(modelF)
modelF$modelsense <- "max"
max_genlate <- gurobi::gurobi(modelF)

bound <- c(min_genlate$objval, max_genlate$objval)

##-------------------------
## Test equivalence of Gamma moments
##-------------------------

test_that("LP problem", {
    expect_equal(result$bound, bound)
    expect_equal(as.numeric(result$lpresult$model$rhs), modelF$rhs)
    expect_equal(result$lpresult$model$sense, modelF$sense)
    expect_equal(as.numeric(result$lpresult$model$A), as.numeric(modelF$A))
    expect_equal(dim(result$lpresult$model$A), dim(modelF$A))
})
