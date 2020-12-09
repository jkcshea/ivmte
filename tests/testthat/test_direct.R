context("Test of direct MTR regression.")
set.seed(10L)

## Only perform this test of Gurobi is available
if (requireNamespace("gurobi", quietly = TRUE)) {
    ## Generate data
    dtm <- ivmte:::gendistMosquito()

    ## -----------------------------------------------------------------------------
    ## Test point identified case
    ## -----------------------------------------------------------------------------

    results <- ivmte(data = dtm,
                     propensity = d ~ 0 + factor(z),
                     m0 = ~ 1 + u + I(u^2),
                     m1 = ~ 1 + u + I(u^2),
                     criterion.tol = 0,
                     target = "ate",
                     treat = 'd',
                     outcome = 'ey',
                     noisy = TRUE)

    ## True coefficients from Mogstad, Torgovitsky (2018, ARE)
    coef.m0 <- c(0.9, -1.1, 0.3)
    coef.m1 <- c(0.35, -0.3, -0.05)
    true.ate <- -0.8/3
    test_that("Point identified case", {
        expect_equal(unname(results$propensity$phat), dtm$pz)
        expect_equal(unname(results$mtr.coef), c(coef.m0, coef.m1))
        expect_equal(results$point.estimate, true.ate)
    })

    ## -----------------------------------------------------------------------------
    ## Test misspecified partially identified case
    ## -----------------------------------------------------------------------------

    resultsAlt <- ivmte(data = dtm,
                        propensity = d ~ 0 + factor(z),
                        m0 = ~ 1 + I(u^2),
                        m1 = ~ 1 + I(u^2),
                        point = FALSE,
                        criterion.tol = 4,
                        target = "ate",
                        treat = 'd',
                        outcome = 'ey',
                        initgrid.nu = 3,
                        audit.nu = 3,
                        m0.inc = TRUE,
                        noisy = TRUE)

    ## Construct design matrix
    dVec <- dtm$d
    pVec <- dtm$pz
    b0.0 <- (monoIntegral(1, 0) - monoIntegral(pVec, 0)) / (1 - pVec)
    b0.1 <- (monoIntegral(1, 1) - monoIntegral(pVec, 1)) / (1 - pVec)
    b0.2 <- (monoIntegral(1, 2) - monoIntegral(pVec, 2)) / (1 - pVec)
    b1.0 <- (monoIntegral(pVec, 0) - monoIntegral(0, 0)) / pVec
    b1.1 <- (monoIntegral(pVec, 1) - monoIntegral(0, 1)) / pVec
    b1.2 <- (monoIntegral(pVec, 2) - monoIntegral(0, 2)) / pVec
    fullY <- dtm$ey
    fullX <- cbind(b0.0 * (1 - dVec),
                   b0.1 * (1 - dVec),
                   b0.2 * (1 - dVec),
                   b1.0 * dVec,
                   b1.1 * dVec,
                   b1.2 * dVec)
    fullFit <- lm.fit(x = fullX, y = fullY)
    ## Check that MTR components are constructed correctly
    test_that("Verify MTR component construction", {
        expect_equal(unname(results$mtr.coef), unname(fullFit$coef))
    })

    ## Now perform regression using misspecified model and ensure that the
    ## coefficients and SSR match with those of the function
    misX <- fullX[, -c(2, 5)]
    misFit <- lm.fit(x = misX, y = fullY)
    misQ <- sum(misFit$resid^2)
    test_that("Verify direct regressions match", {
        expect_equal(unname(resultsAlt$s.set$init.coef), unname(misFit$coef))
        expect_equal(resultsAlt$s.set$Q, misQ)
    })

    ## Now perform the audit using the full grid, as specified in the call
    ## above.
    uGrid <- seq(0, 1, 0.25)
    misQuad <- function(u) {
        c(1, u^2)
    }
    ## Construct base matrices for boundedness
    A0 <- A1 <- t(sapply(uGrid, misQuad))
    A <- rbind(cbind(A0, matrix(0, ncol = 2, nrow = 5)),
               cbind(matrix(0, ncol = 2, nrow = 5), A1))
    ## Duplicate base matrices: one for lb, another for ub
    A <- rbind(A, A)
    ## Construct the monotonicity constraint matrices
    monoA0 <- A0[-1, ] - A0[-5, ]
    monoA <- cbind(monoA0, matrix(0, ncol = 2, nrow = 4))
    ## Combine all constraint matrices together
    A <- rbind(A, monoA)
    ## Define the sense and rhs vectors
    sense <- c(rep('>=', 10), rep('<=', 10), rep('>=', 4))
    rhs <- c(rep(min(fullY), 10), rep(max(fullY), 10), rep(0, 4))
    ## Construct target parameter MTR vector
    tau <- c(-(monoIntegral(1, 0) - monoIntegral(0, 0)),
             -(monoIntegral(1, 2) - monoIntegral(0, 2)),
             (monoIntegral(1, 0) - monoIntegral(0, 0)),
             (monoIntegral(1, 2) - monoIntegral(0, 2)))
    ## Now define quadratic constraints
    yy <- sum(fullY^2)
    q <- -2 * t(misX) %*% fullY
    Q <- t(misX) %*% misX
    ## Now construct the object for gurobi
    model <- list()
    model$obj <- tau
    model$A <- A
    model$sense <- sense
    model$rhs <- rhs
    model$lb <- rep(-Inf, 4)
    model$ub <- rep(Inf, 4)
    qc <- list()
    qc$q <- as.vector(q)
    qc$Qc <- Q
    qc$rhs <- misQ * (1 + 32) - yy
    model$quadcon <- list(qc)
    ## Optimize
    model$modelsense <- 'min'
    qpMin <- gurobi::gurobi(model, list(nonconvex = 0))
    model$modelsense <- 'max'
    qpMax <- gurobi::gurobi(model, list(nonconvex = 0))
    bounds <- c(qpMin$objval, qpMax$objval)
    mtrMin <- qpMin$x
    mtrMax <- qpMax$x

    ## Now check that bounds and MTR coefficients match
    test_that("Verify bounds and MTR coefficient estimates", {
        expect_equal(unname(resultsAlt$bounds), bounds)
        expect_equal(unname(c(resultsAlt$gstar.coef$min.g0,
                              resultsAlt$gstar.coef$min.g1)),
                     mtrMin)
        expect_equal(unname(c(resultsAlt$gstar.coef$max.g0,
                              resultsAlt$gstar.coef$max.g1)),
                     mtrMax)
    })
}
