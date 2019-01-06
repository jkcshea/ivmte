context("Test of case involving splines.")

set.seed(10L)

##------------------------
## Run MST estimator
##------------------------

ivlike <- c(ey ~ d,
            ey ~ d + x,
            ey ~ d + x | z + x)

components <- l(c(intercept, d), d, c(d, x))

set.seed(10L)
result <- ivmte(ivlike = ivlike,
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
                obseq.tol = 0.01,
                grid.nu = 3,
                grid.nx = 2,
                audit.nx = 1,
                audit.nu = 5,
                m1.ub = 55,
                m0.lb = 0,
                mte.inc = TRUE,
                lpsolver = "lpSolveAPI")

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
gstar <- genGammaSplinesTT(distr = dts,
                            u1s1 = t1s1,
                            u0s1 = t0s1,
                            u0s2 = t0s2,
                            weight = wAttSplines,
                            ed = m[["p"]],
                            target = TRUE)
gstar$g0 <- - gstar$g0

## S-set Gamma terms

g1 <- genGammaSplinesTT(distr = dts,
                         u1s1 = u1s1,
                         u0s1 = u0s1,
                         u0s2 = u0s2,
                         weight = sOlsSplines,
                         j = 1,
                         exx = exx[1:2, 1:2])

g2 <- genGammaSplinesTT(distr = dts,
                         u1s1 = u1s1,
                         u0s1 = u0s1,
                         u0s2 = u0s2,
                         weight = sOlsSplines,
                         j = 2,
                         exx = exx[1:2, 1:2])

g3 <- genGammaSplinesTT(distr = dts,
                         u1s1 = u1s1,
                         u0s1 = u0s1,
                         u0s2 = u0s2,
                         zvars = "x",
                         weight = sOlsSplines,
                         j = 2,
                         exx = exx)

g4 <- genGammaSplinesTT(distr = dts,
                         u1s1 = u1s1,
                         u0s1 = u0s1,
                         u0s2 = u0s2,
                         zvars = c("z", "x"),
                         weight = sTslsSplines,
                         j = 2,
                         exz = exz,
                         pi = pi)

g5 <- genGammaSplinesTT(distr = dts,
                         u1s1 = u1s1,
                         u0s1 = u0s1,
                         u0s2 = u0s2,
                         zvars = c("z", "x"),
                         weight = sTslsSplines,
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

## Construct full Gurobi/lpSolveAPI model
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
## modelO$modelsense <- "min"
## minobseq <- gurobi::gurobi(modelO)$objval
minobseq <- runLpSolveAPI(modelO, 'min')$objval

##-------------------------
## Obtain the bounds for the ATT
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

## Code for Gurobi:
## modelF$modelsense <- "min"
## minAtt <- gurobi::gurobi(modelF)
## modelF$modelsense <- "max"
## maxAtt <- gurobi::gurobi(modelF)

minAtt <- runLpSolveAPI(modelF, 'min')
maxAtt <- runLpSolveAPI(modelF, 'max')

bound <- c(minAtt$objval, maxAtt$objval)

##-------------------------
## Test bounds
##-------------------------

test_that("LP problem", {
    expect_equal(result$bound, bound)
    expect_equal(as.numeric(result$lpresult$model$rhs), modelF$rhs)
    expect_equal(result$lpresult$model$sense, modelF$sense)
    expect_equal(as.numeric(result$lpresult$model$A), as.numeric(modelF$A))
    expect_equal(dim(result$lpresult$model$A), dim(modelF$A))
})

##------------------------
## Alternative, where custom spline weights are declared
##------------------------

## Declare weight functions and knots
weight11 <- function(z) {
    1 / (1 + z)
}

weight01 <- function(z, x) {
    - 1 / (1 + z + abs(x))
}
weight02 <- function(z) {
    - (1 / (1 + z)) * 1.33
}

knots01 <- function(z, x) {
    0.3 + 0.3 * z + 0.1 * x
}

## Custom weight estimate using smaller sample
set.seed(10L)
resultAlt <- ivmte(ivlike = ivlike,
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
                   target.weight0 = c(weight01, weight02),
                   target.weight1 = weight11,
                   target.knots0 = knots01,
                   obseq.tol = 0.01,
                   grid.nu = 3,
                   grid.nx = 2,
                   audit.nx = 1,
                   audit.nu = 5,
                   m1.ub = 55,
                   m0.lb = 0,
                   mte.inc = TRUE,
                   lpsolver = "lpSolveAPI")

##------------------------
## Reconstruct target gamma moments
##------------------------

## Construct weights
dts$weight11 <- 1 / (1 + dts$z)
dts$weight01 <- - 1 / (1 + dts$z + abs(dts$x))
dts$weight02 <- - (1 / (1 + dts$z)) * 1.33
dts$knots01 <- 0.3 + 0.3 * dts$z + 0.1 * dts$x

## Construct control group, nonsplines term
u2vecA <- sapply(dts$knots01, function(x) (1 / 3) * x ^ 3)
u2vecB <- sapply(rep(1, nrow(dts)), function(x) (1 / 3) * x ^ 3) - u2vecA

u2vecA <- u2vecA * dts$weight01 * dts$f
u2vecB <- u2vecB * dts$weight02 * dts$f

gstarCustomNS0 <- sum(u2vecA + u2vecB)

## Construct control group, splines 1
splineMat0a <- t(mapply(splines2::ibs,
                        x = dts$knots01,
                        MoreArgs = list(degree = 0,
                                        knots = c(0.2, 0.5, 0.8),
                                        intercept = TRUE,
                                        Boundary.knots = c(0, 1))))

splineMat0b <- t(mapply(splines2::ibs,
                        x = rep(1, nrow(dts)),
                        MoreArgs = list(degree = 0,
                                        knots = c(0.2, 0.5, 0.8),
                                        intercept = TRUE,
                                        Boundary.knots = c(0, 1))))

splineMat01a <- sweep(splineMat0a, MARGIN = 1, STATS = dts$x, FUN = "*")
splineMat01a <- sweep(splineMat01a, MARGIN = 1, STATS = dts$weight01, FUN = "*")
splineMat01a <- sweep(splineMat01a, MARGIN = 1, STATS = dts$f, FUN = "*")

splineMat01b <- splineMat0b - splineMat0a
splineMat01b <- sweep(splineMat01b, MARGIN = 1, STATS = dts$x, FUN = "*")
splineMat01b <- sweep(splineMat01b, MARGIN = 1, STATS = dts$weight02, FUN = "*")
splineMat01b <- sweep(splineMat01b, MARGIN = 1, STATS = dts$f, FUN = "*")

gstarCustomS01 <- colSums(splineMat01a + splineMat01b)

## Construct control group, splines 2

splineMat0c <- t(mapply(splines2::ibs,
                        x = dts$knots01,
                        MoreArgs = list(degree = 1,
                                        knots = c(0.4),
                                        intercept = TRUE,
                                        Boundary.knots = c(0, 1))))

splineMat0d <- t(mapply(splines2::ibs,
                        x = rep(1, nrow(dts)),
                        MoreArgs = list(degree = 1,
                                        knots = c(0.4),
                                        intercept = TRUE,
                                        Boundary.knots = c(0, 1))))

splineMat01c <- sweep(splineMat0c, MARGIN = 1, STATS = dts$weight01, FUN = "*")
splineMat01c <- sweep(splineMat01c, MARGIN = 1, STATS = dts$f, FUN = "*")

splineMat01d <- splineMat0d - splineMat0c
splineMat01d <- sweep(splineMat01d, MARGIN = 1, STATS = dts$weight02, FUN = "*")
splineMat01d <- sweep(splineMat01d, MARGIN = 1, STATS = dts$f, FUN = "*")

gstarCustomS02 <- colSums(splineMat01c + splineMat01d)

## Construct treated  group, nonsplines term

gstarCustomNS1 <- c(sum(dts$weight11 * dts$f),
                    sum(dts$x * dts$weight11 * dts$f))

## Construct treated  group, splines

splineMat1 <- t(mapply(splines2::ibs,
                        x = rep(1, nrow(dts)),
                        MoreArgs = list(degree = 2,
                                        knots = c(0.3, 0.6),
                                        intercept = FALSE,
                                        Boundary.knots = c(0, 1))))

splineMat1 <- sweep(splineMat1, MARGIN = 1, STATS = dts$weight11, FUN = "*")
splineMat1 <- sweep(splineMat1, MARGIN = 1, STATS = dts$f, FUN = "*")

gstarCustomS11 <- colSums(splineMat1)

## Group terms
gstarCustom0 <- c(gstarCustomNS0,
                  gstarCustomS02,
                  gstarCustomS01)

gstarCustom1 <- c(gstarCustomNS1,
                  gstarCustomS11)

## Begin testing
test_that("Custom weights", {
    expect_equal(as.numeric(gstarCustom0), as.numeric(resultAlt$gstar$g0))
    expect_equal(as.numeric(gstarCustom1), as.numeric(resultAlt$gstar$g1))
})
