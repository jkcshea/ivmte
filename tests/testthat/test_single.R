context("Test estimation when only a single IV-like specification is provided.")
set.seed(10L)

##------------------------
## Run MST estimator
##------------------------

dtcf <- ivmte:::gendistCovariates()$data.full
dtc <- ivmte:::gendistCovariates()$data.dist
result <- ivmte(ivlike = ey ~ 1 +d + x1 + x2,
                data = dtcf,
                components = l(d, x1),
                subset = z2 %in% c(2, 3),
                propensity = d ~ x1 + x2 + z1 + z2,
                link = "logit",
                m0 = ~ x1 + x2:u + x2:I(u^2),
                m1 = ~ x1 + x1:x2 + u + x1:u + x2:I(u^2),
                uname = u,
                target = "late",
                late.from = c(z1 = 1, z2 = 2),
                late.to = c(z1 = 0, z2 = 3),
                late.X = c(x1 = 0, x2 = 1),
                criterion.tol = 0.01,
                initgrid.nu = 4,
                initgrid.nx = 2,
                audit.nx = 5,
                audit.nu = 5,
                solver = "lpSolveAPI")

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

##-------------------------
## Construct OLS estimate 3, with controls and subsetting
##-------------------------

mv <- popmean(varlist, subset(dtc, dtc$z2 %in% c(2, 3)))
m <- as.list(mv)
names(m) <- rownames(mv)
exx <- symat(c(1, m[["p"]], m[["x1"]], m[["x2"]],
               m[["p"]], m[["I(p * x1)"]], m[["I(p * x2)"]],
               m[["I(x1 * x1)"]], m[["I(x1 * x2)"]],
               m[["I(x2 * x2)"]]))
exy <- matrix(c(m[["ey"]], m[["eyd"]],
                m[["I(ey * x1)"]],
                m[["I(ey * x2)"]]))
ols <- (solve(exx) %*% exy)

##-------------------------
## Test equivalence of IV-like estimates
##-------------------------

test_that("IV-like estimates", {
    expect_equal(as.numeric(result$s.set$s1$beta), as.numeric(ols[2]))
    expect_equal(as.numeric(result$s.set$s2$beta), as.numeric(ols[3]))
})

##-------------------------
## Construct gamma terms for OLS, with controls and subset
##-------------------------

## Generate weights
dtc.x <- split(as.matrix(dtc[, c("x1", "x2")]), seq(1, nrow(dtc)))
## Replace propensity score with the fitted values
fit <- glm(d ~ x1 + x2 + z1 + z2, family = binomial(link = "logit"),
           data = dtcf)
dtc$p <- predict(fit, dtc, type = "response")
dtc$s.ols.0.d <- unlist(lapply(dtc.x, sOls3, d = 0, j = 2, exx = exx))
dtc$s.ols.1.d <- unlist(lapply(dtc.x, sOls3, d = 1, j = 2, exx = exx))
dtc$s.ols.0.x1 <- unlist(lapply(dtc.x, sOls3, d = 0, j = 3, exx = exx))
dtc$s.ols.1.x1 <- unlist(lapply(dtc.x, sOls3, d = 1, j = 3, exx = exx))
g.ols.d  <- genGammaTT(subset(dtc, dtc$z2 %in% c(2, 3)),
                      "s.ols.0.d",
                      "s.ols.1.d")
g.ols.x1 <- genGammaTT(subset(dtc, dtc$z2 %in% c(2, 3)),
                      "s.ols.0.x1",
                      "s.ols.1.x1")
## NOTE: the term for the constants are approximately 0, 1e-16.

##-------------------------
## Construct target gammas
##-------------------------

## LATE for fixed X
late.ub <- subset(dtc,
                  dtc$z1 == 0 &
                  dtc$z2 == 3 &
                  dtc$x1 == 0 &
                  dtc$x2 == 1)$p
late.lb <- subset(dtc,
                  dtc$z1 == 1 &
                  dtc$z2 == 2 &
                  dtc$x1 == 0 &
                  dtc$x2 == 1)$p
dtc$w.late.1 <- 1 / (late.ub - late.lb)
dtc$w.late.0 <- - dtc$w.late.1
g.star.late <- genGammaTT(dtc[dtc$x1 == 0 & dtc$x2 == 1, ],
                        "w.late.0",
                        "w.late.1",
                        lb = late.lb,
                        ub = late.ub)

##-------------------------
## Test equivalence of Gamma terms
##-------------------------

test_that("Gamma moments", {
    expect_equal(as.numeric(c(result$gstar$g0, result$gstar$g1)),
                 as.numeric(unlist(g.star.late)))
    expect_equal(as.numeric(c(result$s.set$s1$g0, result$s.set$s1$g1)),
                 as.numeric(unlist(g.ols.d)))
    expect_equal(as.numeric(c(result$s.set$s2$g0, result$s.set$s2$g1)),
                 as.numeric(unlist(g.ols.x1)))
})

##-------------------------
## Minimize observational equivalence
##-------------------------

estimates <- c(ols[c(2, 3)])
## Construct A matrix components
A <- rbind(c(g.ols.d$g0, g.ols.d$g1),
           c(g.ols.x1$g0, g.ols.x1$g1))
Aextra <- matrix(0, nrow = nrow(A), ncol = 2 * nrow(A))
for (i in 1:nrow(A)) {
    Aextra[i, (i * 2 - 1)] <- -1
    Aextra[i, (i * 2)] <- 1
}
grid <- result$audit.grid$initial[, 1:3]


xGrid <- result$audit.grid$audit.x
nx <- nrow(xGrid)
uGrid <- result$audit.grid$audit.u
xGrid <- xGrid[rep(seq(nrow(xGrid)), each = length(uGrid)), ]
uGrid <- rep(uGrid, times = nx)
grid <- cbind(xGrid, uGrid)
colnames(grid) <- c("x1", "x2", "u")
grid <- data.frame(grid)
mono0 <- model.matrix(~ x1 + x2:u + x2:I(u^2),
                      data = grid)
mono1 <- model.matrix(~ x1 + x1:x2 + u + x1:u + x2:I(u^2),
                      data = grid)

## Construct boundedness matrix components
maxy <- max(subset(dtc, dtc$z2 %in% c(2, 3))[, c("ey0", "ey1")])
miny <- min(subset(dtc, dtc$z2 %in% c(2, 3))[, c("ey0", "ey1")])
Bzeroes <- matrix(0, ncol = ncol(Aextra), nrow(grid))
b0zeroes <- matrix(0, ncol = ncol(mono0), nrow = nrow(grid))
b1zeroes <- matrix(0, ncol = ncol(mono1), nrow = nrow(grid))
m0bound <- cbind(Bzeroes, mono0, b1zeroes)
m1bound <- cbind(Bzeroes, b0zeroes, mono1)
mtebound <- cbind(Bzeroes, -mono0, mono1)

## Construct full lpSolveAPI model
modelO <- list()
modelO$obj <- c(replicate(ncol(Aextra), 1),
                replicate(ncol(A), 0))
modelO$rhs <- c(estimates,
                replicate(nrow(m0bound), miny),
                replicate(nrow(m1bound), miny),
                replicate(nrow(mtebound), miny - maxy),
                replicate(nrow(m0bound), maxy),
                replicate(nrow(m1bound), maxy),
                replicate(nrow(mtebound), maxy - miny))
modelO$sense <- c(replicate(length(estimates), "="),
                  replicate(nrow(m0bound), ">="),
                  replicate(nrow(m1bound), ">="),
                  replicate(nrow(mtebound), ">="),
                  replicate(nrow(m0bound), "<="),
                  replicate(nrow(m1bound), "<="),
                  replicate(nrow(mtebound), "<="))
modelO$A <- rbind(cbind(Aextra, A),
                  m0bound,
                  m1bound,
                  mtebound,
                  m0bound,
                  m1bound,
                  mtebound)
modelO$ub <- c(replicate(ncol(Aextra), Inf),
                replicate(ncol(A), Inf))
modelO$lb <- c(replicate(ncol(Aextra), 0),
                replicate(ncol(A), -Inf))
## Minimize observational equivalence deviation
lpsolver.options <- list(epslevel = "tight")
minobseq <- runLpSolveAPI(modelO, 'min', lpsolver.options)$objval

##-------------------------
## Obtain the bounds for the LATE
##-------------------------

tolerance <- 1.01
Atop <- c(replicate(ncol(Aextra), 1),
          replicate(ncol(A), 0))
modelF <- list()
modelF$obj <- c(replicate(ncol(Aextra), 0),
                g.star.late$g0,
                g.star.late$g1)
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
minLate <- runLpSolveAPI(modelF, 'min', lpsolver.options)
maxLate <- runLpSolveAPI(modelF, 'max', lpsolver.options)
bound <- c(lower = minLate$objval, upper = maxLate$objval)

##-------------------------
## Test bounds
##-------------------------

test_that("LP problem", {
    expect_equal(result$bound, bound)
})
