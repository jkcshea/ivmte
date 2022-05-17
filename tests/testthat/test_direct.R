context("Test of direct MTR regression.")
set.seed(10L)

## Only perform this test of Gurobi is available
if (requireNamespace("gurobi", quietly = TRUE)) {
    ## Generate data
    dtm <- ivmte:::gendistMosquito()

    ## -------------------------------------------------------------------------
    ## Test point identified case
    ## -------------------------------------------------------------------------

    results <- ivmte(data = dtm,
                     propensity = d ~ 0 + factor(z),
                     m0 = ~ 1 + u + I(u^2),
                     m1 = ~ 1 + u + I(u^2),
                     criterion.tol = 0,
                     target = "ate",
                     treat = 'd',
                     outcome = 'ey',
                     noisy = TRUE,
                     solver = 'gurobi')

    ## True coefficients from Mogstad, Torgovitsky (2018, ARE)
    coef.m0 <- c(0.9, -1.1, 0.3)
    coef.m1 <- c(0.35, -0.3, -0.05)
    true.ate <- -0.8/3
    test_that("Point identified case", {
        expect_equal(unname(results$propensity$phat), dtm$pz)
        expect_equal(unname(results$mtr.coef), c(coef.m0, coef.m1))
        expect_equal(results$point.estimate, true.ate)
    })

    ## -------------------------------------------------------------------------
    ## Test misspecified partially identified case
    ## -------------------------------------------------------------------------

    ## Test misspecified case/partially identified case
    criterion.tol <- 0.5
    dtm$x <- 1
    devtools::load_all("../ivmte")
    resultsAlt <- ivmte(data = dtm,
                        propensity = d ~ 0 + factor(z),
                        m0 = ~ 1 + u + I(u^2),
                        m1 = ~ 1 + u + I(u^2) + x,
                        point = TRUE,
                        criterion.tol = criterion.tol,
                        target = 'ate',
                        outcome = 'ey',
                        initgrid.nu = 3,
                        audit.nu = 3,
                        m0.inc = TRUE,
                        m1.inc = TRUE,
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

    ## Now perform regression using misspecified collinear model.
    misX <- cbind(fullX, dtm$x * dtm$d)
    ## misX <- cbind(fullX, dtm$x * dtm$d)
    tmp.y <- dtm$ey
    eyb <- NULL
    ebb <- NULL
    for (i in 1:ncol(misX)) {
        tmp.b <- misX[, i]
        tmp.bmat <- colMeans(sweep(misX, 1, tmp.b, "*"))
        eyb <- c(eyb, mean(tmp.b * tmp.y))
        ebb <- rbind(ebb, tmp.bmat)
    }
    ## Now generate the shape constraints.
    uGrid <- seq(0, 1, 0.25)
    evalQuad <- function(u) {
        c(1, u, u^2)
    }
    evalColQuad <- function(u) {
        c(1, u, u^2, 1)
    }
    ## Construct base matrices for boundedness
    A0 <- t(sapply(uGrid, evalQuad))
    A1 <- t(sapply(uGrid, evalColQuad))
    A <- rbind(cbind(A0, matrix(0, ncol = 4, nrow = 5)),
               cbind(matrix(0, ncol = 3, nrow = 5), A1))
    ## Duplicate base matrices: one for lb, another for ub
    A <- rbind(A, A)
    ## Construct the monotonicity constraint matrices
    monoA0 <- A0[-1, ] - A0[-5, ]
    monoA0 <- cbind(monoA0, matrix(0, ncol = 4, nrow = 4))
    monoA1 <- A1[-1, ] - A1[-5, ]
    monoA1 <- cbind(matrix(0, ncol = 3, nrow = 4), monoA1)
    ## Combine all constraint matrices together
    A <- rbind(ebb, A, monoA0, monoA1)
    ## Define the sense and rhs vectors
    sense <- c(rep('=', 7), rep('>=', 10), rep('<=', 10), rep('>=', 8))
    rhs <- c(eyb, rep(min(fullY), 10), rep(max(fullY), 10), rep(0, 8))
    ## Now add in the slack variables
    Aextra <- matrix(0, nrow = nrow(A), ncol = 2 * ncol(misX))
    for (i in 1:ncol(misX)) {
        Aextra[i, (i * 2 - 1)] <- -1
        Aextra[i, (i * 2)] <- 1
    }
    ## Add lower and upper bound resrictions
    ub <- c(replicate(ncol(Aextra), Inf),
            replicate(ncol(A), Inf))
    lb <- c(replicate(ncol(Aextra), 0),
            replicate(ncol(A), -Inf))
    ## Define objective
    obj <- c(rep(1, ncol(Aextra)), rep(0, ncol(A)))
    ## Now combine everything together
    A <- cbind(Aextra, A)
    ## Now combine the model together
    model.criterion <- list(obj = obj,
                            A = A,
                            rhs = rhs,
                            sense = sense,
                            lb = lb,
                            ub = ub,
                            modelsense = "min")
    ## Check that minimum criterion matches
    minobseq <- gurobi::gurobi(model.criterion)
    test_that("Verify criterions match", {
        expect_equal(minobseq$objval, resultsAlt$audit.criterion)
    })
    ## Note: Note that even if the model is correctly specified, the
    ## minimum criterion would not be zero in this example. This is
    ## because of the shape constraints. The data set 'dtm' is
    ## generated using mean counterfactual outcomes, so it will be
    ## possible to match all the moment restrictions. However, the
    ## shape restrictions are not based on moments, but on specific
    ## values of the unobserved variable 'u'.

    ## Update model with minimum criterion
    A.criterion <- c(rep(1, ncol(Aextra)), rep(0, ncol(A0) + ncol(A1)))
    rhs.criterion <- minobseq$objval * (1 + criterion.tol)
    sense.criterion <- "<"
    A <- rbind(A.criterion, A)
    rhs <- c(rhs.criterion, rhs)
    sense <- c(sense.criterion, sense)
    ## Update objective
    t0.0 <- -monoIntegral(1, 0)
    t0.1 <- -monoIntegral(1, 1)
    t0.2 <- -monoIntegral(1, 2)
    t1.0 <- monoIntegral(1, 0)
    t1.1 <- monoIntegral(1, 1)
    t1.2 <- monoIntegral(1, 2)
    t1.x <- mean(dtm$x)
    obj.te <- c(t0.0, t0.1, t0.2, t1.0, t1.1, t1.2, t1.x)
    model <- copy(model.criterion)
    model$obj <- c(rep(0, ncol(Aextra)), obj.te)
    model$A <- A
    model$rhs <- rhs
    model$sense <- sense
    ## Obtain lower bounds and check for equality
    result.min <- gurobi::gurobi(model)
    model$modelsense <- "max"
    result.max <- gurobi::gurobi(model)
    bounds <- c(result.min$objval, result.max$objval)
    test_that("Verify bounds and MTR coefficient estimates", {
        expect_equal(unname(resultsAlt$bounds), bounds)
    })
}
