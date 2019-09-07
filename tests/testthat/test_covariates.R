context("Test of case involving only covariates, no splines.")
set.seed(10L)

##------------------------
## Run MST estimator
##------------------------

dtcf <- ivmte:::gendistCovariates()$data.full
ivlike <- c(ey ~ d,
            ey ~ d + x1,
            ey ~ d + x1 + x2,
            ey ~ d + x1 + x2 | x1 + x2 + z1 + z2,
            ey ~ d | factor(z2))
components <- l(d, d, c(d, x1, x2), d, d)
subsets    <- l(, , z2 %in% c(2, 3), , z2 %in% c(2, 3))
result <- ivmte(ivlike = ivlike,
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
                obseq.tol = 0.01,
                initgrid.nu = 3,
                initgrid.nx = 2,
                audit.nx = 1,
                audit.nu = 5,
                m0.inc = TRUE,
                m1.inc = TRUE,
                mte.dec = TRUE,
                treat = d,
                lpsolver = "lpSolveAPI",
                seed = 10L)

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

ey.z2.2 <- popmean(~ 0 + ey, subset(dtc, dtc$z2 == 2))
ey.z2.3 <- popmean(~ 0 + ey, subset(dtc, dtc$z2 == 3))
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
dtc$s.ols1.0.d <- sOls1d(0, exx = ols1.exx)
dtc$s.ols1.1.d <- sOls1d(1, exx = ols1.exx)

## Gammas for D = 0
## m0 = ~ x1 + I(x2 * u) + I(x2 * u^2),
## ols1.0.d.0 means "OLS specification 2. For D = 0. For variable
## "d". For term 0 in md."
g.ols1 <- genGammaTT(dtc, "s.ols1.0.d", "s.ols1.1.d")

##-------------------------
## Construct gamma terms for OLS, with single control
##-------------------------

## Generate weights
dtc$s.ols2.0.d <- sapply(dtc$x1, sOls2d,
                         d = 0,
                         exx = ols1.exx)
dtc$s.ols2.1.d <- sapply(dtc$x1, sOls2d,
                         d = 1,
                         exx = ols1.exx)
g.ols2 <- genGammaTT(dtc, "s.ols2.0.d", "s.ols2.1.d")

##-------------------------
## Construct gamma terms for OLS, with controls and subset
##-------------------------

## Generate weights
dtc.x <- split(as.matrix(dtc[, c("x1", "x2")]), seq(1, nrow(dtc)))
dtc$s.ols3.0.d <- unlist(lapply(dtc.x, sOls3,
                                d = 0,
                                j = 2,
                                exx = ols2.exx))
dtc$s.ols3.1.d <- unlist(lapply(dtc.x, sOls3,
                                d = 1,
                                j = 2,
                                exx = ols2.exx))
dtc$s.ols3.0.x1 <- unlist(lapply(dtc.x, sOls3,
                                 d = 0,
                                 j = 3,
                                 exx = ols2.exx))
dtc$s.ols3.1.x1 <- unlist(lapply(dtc.x, sOls3,
                                 d = 1,
                                 j = 3,
                                 exx = ols2.exx))
dtc$s.ols3.0.x2 <- unlist(lapply(dtc.x, sOls3,
                                 d = 0,
                                 j = 4,
                                 exx = ols2.exx))
dtc$s.ols3.1.x2 <- unlist(lapply(dtc.x, sOls3,
                                 d = 1,
                                 j = 4,
                                 exx = ols2.exx))
g.ols3.d  <- genGammaTT(subset(dtc, dtc$z2 %in% c(2, 3)),
                        "s.ols3.0.d",
                        "s.ols3.1.d")
g.ols3.x1 <- genGammaTT(subset(dtc, dtc$z2 %in% c(2, 3)),
                        "s.ols3.0.x1",
                        "s.ols3.1.x1")
g.ols3.x2 <- genGammaTT(subset(dtc, dtc$z2 %in% c(2, 3)),
                        "s.ols3.0.x2",
                        "s.ols3.1.x2")
## NOTE: the term for the constants are approximately 0, 1e-16.

##-------------------------
## Construct gamma terms for TSLS
##-------------------------

## Generate weights
dtc.z <- split(as.matrix(dtc[, c("z1", "z2", "x1", "x2")]), seq(1, nrow(dtc)))
dtc$s.tsls.0.d <- unlist(lapply(dtc.z, sTsls,
                                j = 2,
                                exz = tsls.exz,
                                pi  = tsls.pi))
dtc$s.tsls.1.d <- unlist(lapply(dtc.z, sTsls,
                                j = 2,
                                exz = tsls.exz,
                                pi  = tsls.pi))
g.tsls <- genGammaTT(dtc, "s.tsls.0.d", "s.tsls.1.d")

##-------------------------
## Construct gamma terms for simple Wald (i.e. one instrument)
##-------------------------

p.z2.2 <- sum(subset(dtc, dtc$z2 == 2)$f)
p.z2.3 <- sum(subset(dtc, dtc$z2 == 3)$f)
## Generate weights
dtc$s.wald.0.d <- sapply(dtc$z2, sWald,
                         p.to   = p.z2.3,
                         p.from = p.z2.2,
                         e.to   = ed.z2.3,
                         e.from = ed.z2.2)
dtc$s.wald.1.d <- sapply(dtc$z2, sWald,
                         p.to   = p.z2.3,
                         p.from = p.z2.2,
                         e.to   = ed.z2.3,
                         e.from = ed.z2.2)
g.wald <- genGammaTT(dtc, "s.wald.0.d", "s.wald.1.d")

##-------------------------
## Construct target gammas
##-------------------------

## Generalized LATE
wald.ub <- 0.7
wald.lb <- 0.2
dtc$w.genlate.1 <- 1 / (wald.ub - wald.lb)
dtc$w.genlate.0 <- - dtc$w.genlate.1
g.star.genlate <- genGammaTT(dtc,
                             "w.genlate.0",
                             "w.genlate.1",
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
grid <- matrix(c(0, 3, 0, 2), nrow = 2, byrow = TRUE)
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

mono0 <- model.matrix(~ x1 + I(x2 * u) + I(x2 * u^2),
                      data = grid)
mono1 <- model.matrix(~ x1 + I(x1 * x2) + u + I(x1 * u) + I(x2 * u^2),
                      data = grid)
monoA0 <- mono0[c(2, 3, 5, 6), ] - mono0[c(1, 2, 4, 5), ]
monoA1 <- mono1[c(2, 3, 5, 6), ] - mono1[c(1, 2, 4, 5), ]
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

##-------------------------
## Expand the above according to the audit
##-------------------------

agrid <- result$audit.grid$audit
amono0 <- model.matrix(~ x1 + I(x2 * u) + I(x2 * u^2),
                      data = agrid)
amono1 <- model.matrix(~ x1 + I(x1 * x2) + u + I(x1 * u) + I(x2 * u^2),
                       data = agrid)
amonoA0 <- amono0[seq(2, 7), ] - amono0[seq(1, 6), ]
amonoA1 <- amono1[seq(2, 7), ] - amono1[seq(1, 6), ]
aAzeroes <- matrix(0, ncol = 14, nrow = 6)
am0zeroes <- matrix(0, ncol = 6, nrow = 6)
am1zeroes <- matrix(0, ncol = 4, nrow = 6)
am0mono  <- cbind(aAzeroes, amonoA0, am0zeroes)
am1mono  <- cbind(aAzeroes, am1zeroes, amonoA1)
amtemono <- cbind(aAzeroes, -amonoA0, amonoA1)

## Construct boundedness matrix components
aBzeroes <- matrix(0, ncol = 14, nrow = nrow(agrid))
ab0zeroes <- matrix(0, ncol = ncol(amono0), nrow = nrow(agrid))
ab1zeroes <- matrix(0, ncol = ncol(amono1), nrow = nrow(agrid))
am0bound <- cbind(aBzeroes, amono0, ab1zeroes)
am1bound <- cbind(aBzeroes, ab0zeroes, amono1)

## Construct the audit matrices
arhs <- c(replicate(nrow(am0bound), miny),
          replicate(nrow(am1bound), miny),
          replicate(nrow(am0bound), maxy),
          replicate(nrow(am1bound), maxy),
          replicate(nrow(am0mono), 0),
          replicate(nrow(am1mono), 0),
          replicate(nrow(amtemono), 0))
asense <- c(replicate(nrow(am0bound), ">="),
            replicate(nrow(am1bound), ">="),
            replicate(nrow(am0bound), "<="),
            replicate(nrow(am1bound), "<="),
            replicate(nrow(am0mono), ">="),
            replicate(nrow(am1mono), ">="),
            replicate(nrow(amtemono), "<="))
aA <- rbind(am0bound,
            am1bound,
            am0bound,
            am1bound,
            am0mono,
            am1mono,
            amtemono)

violateVec <- c(14, 22, 34, 40, 46, 29, 42)
addShapeRhs <- arhs[violateVec]
addShapeSense <- asense[violateVec]
addShapeA <- aA[violateVec, ]

##-------------------------
## Obtain minimum criteiron
##-------------------------

## Construct full lpSolveAPI model
model.o <- list()
model.o$obj <- c(replicate(14, 1), replicate(10, 0))
model.o$rhs <- c(estimates,
                 replicate(nrow(m0bound), miny),
                 replicate(nrow(m1bound), miny),
                 replicate(nrow(m0bound), maxy),
                 replicate(nrow(m1bound), maxy),
                 replicate(nrow(m0mono), 0),
                 replicate(nrow(m1mono), 0),
                 replicate(nrow(mtemono), 0),
                 addShapeRhs)
model.o$sense <- c(replicate(7, "="),
                   replicate(nrow(m0bound), ">="),
                   replicate(nrow(m1bound), ">="),
                   replicate(nrow(m0bound), "<="),
                   replicate(nrow(m1bound), "<="),
                   replicate(nrow(m0mono), ">="),
                   replicate(nrow(m1mono), ">="),
                   replicate(nrow(mtemono), "<="),
                   addShapeSense)
model.o$A <- rbind(cbind(A.extra, A),
                   m0bound,
                   m1bound,
                   m0bound,
                   m1bound,
                   m0mono,
                   m1mono,
                   mtemono,
                   addShapeA)
model.o$ub <- c(replicate(14, Inf), replicate(10, Inf))
model.o$lb <- c(replicate(14, 0), replicate(10, -Inf))
## Minimize observational equivalence deviation
minobseq <- runLpSolveAPI(model.o, 'min')$objval

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
min_genlate <- runLpSolveAPI(model.f, 'min')
max_genlate <- runLpSolveAPI(model.f, 'max')
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
