#' Constructing LP problem
#'
#' This function takes in the IV estimates from the set of IV
#' regressions declared by the user, as well as their corresponding
#' moments of the terms in the MTR. These are then used to construct
#' the components that make up the LP problem. Additional constraint
#' matrix is added using \code{mbA} (\code{mb} stands for
#' "monotonicity/boundedness"); extra model sense is added using
#' \code{mbs}; extra RHS values added using \code{mbrhs}). Depending
#' on the linear programming solver used, this function will return
#' different output specific to the solver.
#' @param sset List of IV-like estimates and the corresponding gamma
#'     terms.
#' @param orig.sset list, only used for bootstraps. The list
#'     caontains the gamma moments for each element in the S-set, as
#'     well as the IV-like coefficients.
#' @param mbA Matrix used to define the constraints in the LP problem.
#' @param mbs Vector of model sense/inequalities signs used to define
#'     the constraints in the LP problem.
#' @param mbrhs Vector of constants used to define the constraints in
#'     the LP problem.
#' @param lpsolver string, name of the package used to solve the LP
#'     problem.
#' @param shape boolean, default set to TRUE. Switch to determine
#'     whether or not to include shape restrictions in the LP problem.
#' @return A list of matrices and vectors necessary to define an LP
#'     problem for Gurobi.
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#'
#' ## Declare empty list to be updated (in the event multiple IV like
#' ## specifications are provided
#' sSet <- list()
#'
#' ## Declare MTR formulas
#' formula1 = ~ 1 + u
#' formula0 = ~ 1 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#' polynomials1 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' ## Generate propensity score model
#' propensityObj <- propensity(formula = d ~ z,
#'                             data = dtm,
#'                             link = "linear")
#'
#' ## Generate target gamma moments
#' ivEstimates <- ivEstimate(formula = ey ~ d | z,
#'                           data = dtm,
#'                           components = l(intercept, d),
#'                           treat = d,
#'                           list = FALSE)
#'
#' ## Construct S-set, which contains the coefficients and weights
#' ## corresponding to various IV-like estimands
#' sSet <- genSSet(data = dtm,
#'                 sset = sSet,
#'                 sest = ivEstimates,
#'                 splinesobj = splinesList,
#'                 pmodobj = propensityObj$phat,
#'                 pm0 = polynomials0,
#'                 pm1 = polynomials1,
#'                 ncomponents = 2,
#'                 scount = 1,
#'                 yvar = "ey",
#'                 dvar = "d",
#'                 means = TRUE)
#'
#' ## Construct the LP problem to be solved using lpSolveAPI
#' lpSetup(sset = sSet$sset, lpsolver = "lpSolveAPI")
#'
#' @export
lpSetup <- function(sset, orig.sset = NULL, mbA = NULL, mbs = NULL,
                    mbrhs = NULL, lpsolver, shape = TRUE) {
    lpsolver <- tolower(lpsolver)
    ## determine lengths
    sn  <- length(sset)
    gn0 <- length(sset$s1$g0)
    gn1 <- length(sset$s1$g1)
    ## generate all vectors/matrices for LP optimization to minimize
    ## observational equivalence
    obj <- c(replicate(sn * 2, 1),
             replicate(gn0 + gn1, 0))
    rhs <- unlist(lapply(sset, function(x) x[["beta"]]))
    if (!is.null(orig.sset)) {
        ## Recenter RHS when bootstrapping
        rhs <- rhs - unlist(lapply(orig.sset, function(x) x[["beta"]]))
    }
    sense <- replicate(sn, "=")
    A <- NULL
    scount <- 0
    for (s in names(sset)) {
        avec <- replicate(2 * sn, 0)
        avec[(2 * scount + 1):(2 * scount + 2)] <- c(-1, 1)
        ## Regarding c(-1, 1), the -1 is for w+, 1 is for w-
        g0fill <- sset[[s]]$g0
        g1fill <- sset[[s]]$g1
        if (!is.null(orig.sset)) {
            ## Recenter gamma vectors when bootstrapping
            g0fill <- g0fill - orig.sset[[s]]$g0
            g1fill <- g1fill - orig.sset[[s]]$g1
        }
        avec <- c(avec, g0fill, g1fill)
        A <- rbind(A, avec)
        scount <- scount + 1
    }
    colnames(A) <- c(seq(1, 2 * sn),
                     colnames(A)[(2 * sn + 1) : ncol(A)])
    ## Define bounds on parameters
    ub <- replicate(ncol(A), Inf)
    lb <- c(replicate(sn * 2, 0), replicate(gn0 + gn1, -Inf))
    ## Add in additional constraints if included
    if (shape == TRUE) {
        mbA     <- rbind(A, mbA)
        sense <- c(sense, mbs)
        rhs   <- c(rhs, mbrhs)
    } else {
        mbA <- A
    }
    rm(A)
    if (lpsolver %in% c("gurobi", "lpsolveapi")) {
        mbA <- Matrix::Matrix(mbA, sparse = TRUE)
    }
    return(list(obj = obj,
                rhs = rhs,
                sense = sense,
                A = mbA,
                ub = ub,
                lb = lb,
                sn = sn,
                gn0 = gn0,
                gn1 = gn1))
}


lpSetupAlt <- function(env, sset, orig.sset = NULL, mbA = NULL, mbs = NULL,
                    mbrhs = NULL, lpsolver, shape = TRUE) {
    lpsolver <- tolower(lpsolver)
    ## determine lengths
    sn  <- length(sset)
    gn0 <- length(sset$s1$g0)
    gn1 <- length(sset$s1$g1)
    ## generate all vectors/matrices for LP optimization to minimize
    ## observational equivalence
    rhs <- unlist(lapply(sset, function(x) x[["beta"]]))
    if (!is.null(orig.sset)) {
        ## Recenter RHS when bootstrapping
        rhs <- rhs - unlist(lapply(orig.sset, function(x) x[["beta"]]))
    }
    sense <- replicate(sn, "=")
    A <- NULL
    scount <- 0
    for (s in names(sset)) {
        avec <- replicate(2 * sn, 0)
        avec[(2 * scount + 1):(2 * scount + 2)] <- c(-1, 1)
        ## Regarding c(-1, 1), the -1 is for w+, 1 is for w-
        g0fill <- sset[[s]]$g0
        g1fill <- sset[[s]]$g1
        if (!is.null(orig.sset)) {
            ## Recenter gamma vectors when bootstrapping
            g0fill <- g0fill - orig.sset[[s]]$g0
            g1fill <- g1fill - orig.sset[[s]]$g1
        }
        avec <- c(avec, g0fill, g1fill)
        A <- rbind(A, avec)
        scount <- scount + 1
    }
    colnames(A) <- c(seq(1, 2 * sn),
                     colnames(A)[(2 * sn + 1) : ncol(A)])
    ## Define bounds on parameters
    ub <- replicate(ncol(A), Inf)
    lb <- c(replicate(sn * 2, 0), replicate(gn0 + gn1, -Inf))
    ## Add in additional constraints if included
    if (shape == TRUE) {
        mbA     <- rbind(A, mbA)
        sense <- c(sense, mbs)
        rhs   <- c(rhs, mbrhs)
    } else {
        mbA <- A
    }
    rm(A)
    if (lpsolver %in% c("gurobi", "lpsolveapi")) {
        mbA <- Matrix::Matrix(mbA, sparse = TRUE)
    }
    env$lpobj <- list(rhs = rhs,
                      sense = sense,
                      A = mbA,
                      ub = ub,
                      lb = lb,
                      sn = sn,
                      gn0 = gn0,
                      gn1 = gn1)
}

lpSetupUnbounded <- function(env, sset) {
    sn <- length(sset)
    env$lpobj$A[1:(2 * sn), ]
    env$lpobj$rhs[1:(2 * sn)]
    env$lpobj$sense[1:(2 * sn)]
}

lpSetupCriterion <- function(env, sset) {
    ## determine lengths
    sn  <- length(sset)
    gn0 <- length(sset$s1$g0)
    gn1 <- length(sset$s1$g1)
    ## generate all vectors/matrices for LP optimization to minimize
    ## observational equivalence
    env$lpobj$obj <- c(replicate(sn * 2, 1),
                       replicate(gn0 + gn1, 0))
}

## This function simply adjust the lpObj so that it is compatible with
## the solver selected.
lpSetupSolver <- function(env, lpsolver) {
    if (lpsolver == "cplexapi") {
        env$lpobj$sense[env$lpobj$sense == "<"]  <- "L"
        env$lpobj$sense[env$lpobj$sense == "<="] <- "L"
        env$lpobj$sense[env$lpobj$sense == ">"]  <- "G"
        env$lpobj$sense[env$lpobj$sense == ">="] <- "G"
        env$lpobj$sense[env$lpobj$sense == "="]  <- "E"
        env$lpobj$sense[env$lpobj$sense == "=="] <- "E"
        env$lpobj$ub[env$lpobj$ub == Inf] <- cplexAPI::CPX_INFBOUND
        env$lpobj$lb[env$lpobj$lb == -Inf] <- -cplexAPI::CPX_INFBOUND
    }
    if (lpsolver == "lpsolveapi") {
        env$lpobj$sense[env$lpobj$sense == "<"]  <- "<="
        env$lpobj$sense[env$lpobj$sense == ">"]  <- ">="
        env$lpobj$sense[env$lpobj$sense == "=="] <- "="
    }
}

lpSetupCriterionBoot <- function(env, sset, orig.sset,
                                 orig.criterion,
                                 criterion.tol = 0, setup = TRUE) {
    if (setup) {
        sn  <- length(sset)
        gn0 <- length(sset$s1$g0)
        gn1 <- length(sset$s1$g1)
        env$lpobj$obj <- c(replicate(sn * 2, 1),
                           replicate(gn0 + gn1, 0))
        ## Prepare to obtain 'recentered' bootstrap
        ## criterion. Specifically, the |S| equality constraints are
        ## centered. Then, the original |S| equality constraints are
        ## added. In addition, 2 * |S| residual variables are added to
        ## the problem. These new residual variables correspond to the
        ## |S| equality constraints from the original, uncentered
        ## sample.
        tmpA <- NULL
        tmpRhs <- NULL
        tmpSense <- NULL
        scount <- 0
        for (s in names(orig.sset)) {
            avec <- replicate(2 * 2 * length(orig.sset), 0)
            avec[(2 * scount + 1):(2 * scount + 2)] <- c(-1, 1)
            avec <- c(avec, orig.sset[[s]]$g0, orig.sset[[s]]$g1)
            tmpA <- rbind(tmpA, avec)
            tmpRhs <- c(tmpRhs, orig.sset[[s]]$beta)
            tmpSense <- c(tmpSense, "=")
            scount <- scount + 1
        }
        avec <- c(rep(1, 2 * length(orig.sset)),
                  rep(0, 2 * length(orig.sset)),
                  rep(0, length(sset$s1$g0) + length(sset$s1$g1)))
        ## Update lpobj
        env$lpobj$ub <- c(rep(Inf, 2 * length(orig.sset)), env$lpobj$ub)
        env$lpobj$lb <- c(rep(0, 2 * length(orig.sset)), env$lpobj$lb)
        env$lpobj$A <- list(a = avec,
                            b = tmpA,
                            c = cbind(matrix(0,
                                             nrow = nrow(env$lpobj$A),
                                             ncol = length(orig.sset) * 2),
                                      env$lpobj$A))
        rm(avec, tmpA)
        env$lpobj$A <- Reduce(rbind, env$lpobj$A)
        env$lpobj$rhs <- c(orig.criterion * (1 + criterion.tol),
                           tmpRhs, env$lpobj$rhs)
        env$lpobj$sense <- c("<=", tmpSense, env$lpobj$sense)
        env$lpobj$obj <- c(rep(0, 2 * length(orig.sset)), env$lpobj$obj)
    } else {
        ## Simply undo the procedure done above.
        removeCol <- 2 * length(orig.sset)
        removeRow <- length(orig.sset) + 1
        env$lpobj$ub <- env$lpobj$ub[-(1:removeCol)]
        env$lpobj$lb <- env$lpobj$lb[-(1:removeCol)]
        env$lpobj$A <- env$lpobj$A[-(1:removeRow), -(1:removeCol)]
        env$lpobj$rhs <- env$lpobj$rhs[-(1:removeRow)]
        env$lpobj$sense <- env$lpobj$sense[-(1:removeRow)]
        env$lpobj$obj <- env$lpobj$obj[-(1:removeCol)]
    }
}

obsEqMinAlt <- function(env, sset, lpsolver, lpsolver.options, debug = FALSE) {
    lpsolver <- tolower(lpsolver)
    if (lpsolver == "gurobi") {
        if (debug & lpsolver.options$outputflag == 1) {
            message("\nMinimum criterion optimization statistics:")
            message("------------------------------------------")
        }
        env$lpobj$modelsense <- "min"
        if (debug) {
            gurobi::gurobi_write(env$lpobj, "lpCriterion.mps")
            model <- env$lpobj
            save(model, file = "lpCriterion.Rdata")
            rm(model)
        }
        result   <- runGurobi(env$lpobj, lpsolver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
        if (debug) cat("\n")
    }
    if (lpsolver == "cplexapi") {
        result   <- runCplexAPI(env$lpobj, cplexAPI::CPX_MIN, lpsolver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
    }
    if (lpsolver == "lpsolveapi") {
        result   <- runLpSolveAPI(env$lpobj, 'min', lpsolver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
    }
    ## provide nicer output
    g0sol <- optx[(2 * env$lpobj$sn + 1) : (2 * env$lpobj$sn + env$lpobj$gn0)]
    g1sol <- optx[(2 * env$lpobj$sn + env$lpobj$gn0 + 1) :
                  (2 * env$lpobj$sn + env$lpobj$gn0 + env$lpobj$gn1)]
    names(g0sol) <- names(sset$gstar$g0)
    names(g1sol) <- names(sset$gstar$g1)
    return(list(obj = obseqmin,
                x = optx,
                g0 = g0sol,
                g1 = g1sol,
                status = status))
    ## Object 'result' will not be returned---unnecessary, and very
    ## memory intensive.
}


lpSetupBound <- function(env, g0, g1, sset, obseq.factor, lpsolver,
                         setup = TRUE) {
    if (setup) {
        lpsolver <- tolower(lpsolver)
        env$lpobj$obj <- c(replicate(2 * env$lpobj$sn, 0), g0, g1)
        env$lpobj$rhs <- c(obseq.factor, env$lpobj$rhs)
        avec <- c(replicate(2 * env$lpobj$sn, 1),
                  replicate(env$lpobj$gn0 + env$lpobj$gn1, 0))
        env$lpobj$A <- rbind(avec, env$lpobj$A)
        if (lpsolver %in% c("gurobi", "lpsolveapi")) {
            env$lpobj$sense <- c("<=", env$lpobj$sense)
        }
        if (lpsolver == "cplexapi") {
            env$lpobj$sense <- c("L", env$lpobj$sense)
        }
    } else {
        env$lpobj$rhs <- env$lpobj$rhs[-1]
        env$lpobj$sense <- env$lpobj$sense[-1]
        env$lpobj$A <- env$lpobj$A[-1, ]
    }
}

boundAlt <- function(env, sset, obseq.factor, lpsolver,
                     lpsolver.options, noisy = FALSE,
                     debug = FALSE) {
    lpsolver <- tolower(lpsolver)
    ## Obtain lower and upper bounds
    if (lpsolver == "gurobi") {
        if (debug & lpsolver.options$outputflag == 1) {
            message("\nLower bound optimization statistics:")
            message("------------------------------------")
        }
        if (debug == TRUE){
            env$lpobj$modelsense <- NULL
            gurobi::gurobi_write(env$lpobj, "lpBound.mps")
            model <- env$lpobj
            save(model, file = "lpBound.Rdata")
            rm(model)
        }
        env$lpobj$modelsense <- "min"
        minresult <- runGurobi(env$lpobj, lpsolver.options)
        min <- minresult$objval
        minstatus <- minresult$status
        minoptx <- minresult$optx
        if (debug & lpsolver.options$outputflag == 1) {
            message("\nUpper bound optimization statistics:")
            message("------------------------------------")
        }
        env$lpobj$modelsense <- "max"
        maxresult <- runGurobi(env$lpobj, lpsolver.options)
        max <- maxresult$objval
        maxstatus <- maxresult$status
        maxoptx <- maxresult$optx
        if (debug) cat("\n")
    }
    if (lpsolver == "cplexapi") {
        minresult <- runCplexAPI(env$lpobj, cplexAPI::CPX_MIN, lpsolver.options)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runCplexAPI(env$lpobj, cplexAPI::CPX_MAX, lpsolver.options)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
    }
    if (lpsolver == "lpsolveapi") {
        minresult <- runLpSolveAPI(env$lpobj, 'min', lpsolver.options)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runLpSolveAPI(env$lpobj, 'max', lpsolver.options)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
    }
    env$lpobj$modelsense <- NULL
    if (maxstatus == 0 || minstatus == 0) {
        return(NULL)
    }
    ming0 <- minoptx[(2 * env$lpobj$sn + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0)]
    ming1 <- minoptx[(2 * env$lpobj$sn + env$lpobj$gn0 + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0 + env$lpobj$gn1)]
    maxg0 <- maxoptx[(2 * env$lpobj$sn + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0)]
    maxg1 <- maxoptx[(2 * env$lpobj$sn + env$lpobj$gn0 + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0 + env$lpobj$gn1)]
    if (hasArg(sset)) {
        names(ming0) <- names(sset$s1$g0)
        names(ming1) <- names(sset$s1$g1)
        names(maxg0) <- names(sset$s1$g0)
        names(maxg1) <- names(sset$s1$g1)
    }
    ## if (noisy) {
    ##     cat("Min status: ", minstatus, "\n", sep = "")
    ##     cat("Max status: ", maxstatus, "\n", sep = "")
    ##     cat("Bound: (", min, ", ", max, ")\n", sep = "")
    ## }
    return(list(max = max,
                maxg0 = maxg0,
                maxg1 = maxg1,
                maxresult = maxresult,
                maxstatus = maxstatus,
                min = min,
                ming0 = ming0,
                ming1 = ming1,
                minresult = minresult,
                minstatus = minstatus))
}

#' Minimizing violation of observational equivalence
#'
#' Given a set of IV-like estimates and the set of matrices/vectors
#' defining an LP problem, this function minimizes the violation of
#' observational equivalence under the L1 norm.
#' @param sset A list of IV-like estimates and the corresponding gamma
#'     terms.
#' @param orig.sset list, only used for bootstraps. The list
#'     caontains the gamma moments for each element in the S-set, as
#'     well as the IV-like coefficients.
#' @param orig.criterion numeric, only used for bootstraps. The scalar
#'     corresponds to the minimum observational equivalence criterion
#'     from the original sample.
#' @param criterion.tol tolerance for violation of observational
#'     equivalence, set to 0 by default.
#' @param lpobj A list of matrices and vectors defining an LP problem.
#' @param lpsolver string, name of the package used to solve the LP
#'     problem.
#' @param lpsolver.options list, each item of the list should
#'     correspond to an option specific to the LP solver
#'     selected.
#' @param debug boolean, indicates whether or not the function should
#'     provide output when obtaining bounds. The option is only
#'     applied when \code{lpsolver = 'gurobi'}. The output provided is
#'     the same as what the Gurobi API would send to the console.
#' @return A list including the minimum violation of observational
#'     equivalence, the solution to the LP problem, and the status of
#'     the solution.
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#'
#' ## Declare empty list to be updated (in the event multiple IV like
#' ## specifications are provided
#' sSet <- list()
#'
#' ## Declare MTR formulas
#' formula1 = ~ 1 + u
#' formula0 = ~ 1 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#' polynomials1 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' ## Generate propensity score model
#' propensityObj <- propensity(formula = d ~ z,
#'                             data = dtm,
#'                             link = "linear")
#'
#' ## Generate IV estimates
#' ivEstimates <- ivEstimate(formula = ey ~ d | z,
#'                           data = dtm,
#'                           components = l(intercept, d),
#'                           treat = d,
#'                           list = FALSE)
#'
#' ## Generate target gamma moments
#' targetGamma <- genTarget(treat = "d",
#'                          m0 = ~ 1 + u,
#'                          m1 = ~ 1 + u,
#'                          target = "atu",
#'                          data = dtm,
#'                          splinesobj = splinesList,
#'                          pmodobj = propensityObj,
#'                          pm0 = polynomials0,
#'                          pm1 = polynomials1,
#'                          point = FALSE)
#'
#' ## Construct S-set. which contains the coefficients and weights
#' ## corresponding to various IV-like estimands
#' sSet <- genSSet(data = dtm,
#'                 sset = sSet,
#'                 sest = ivEstimates,
#'                 splinesobj = splinesList,
#'                 pmodobj = propensityObj$phat,
#'                 pm0 = polynomials0,
#'                 pm1 = polynomials1,
#'                 ncomponents = 2,
#'                 scount = 1,
#'                 yvar = "ey",
#'                 dvar = "d",
#'                 means = TRUE)
#'
#' ## Define additional upper- and lower-bound constraints for the LP
#' ## problem
#' A <- matrix(0, nrow = 22, ncol = 4)
#' A <- cbind(A, rbind(cbind(1, seq(0, 1, 0.1)),
#'                     matrix(0, nrow = 11, ncol = 2)))
#' A <- cbind(A, rbind(matrix(0, nrow = 11, ncol = 2),
#'                     cbind(1, seq(0, 1, 0.1))))
#'
#' sense <- c(rep(">", 11), rep("<", 11))
#' rhs <- c(rep(0.2, 11), rep(0.8, 11))
#'
#' ## Construct LP object to be interpreted and solved by lpSolveAPI
#' lpObject <- lpSetup(sset = sSet$sset,
#'                     mbA = A,
#'                     mbs = sense,
#'                     mbrhs = rhs,
#'                     lpsolver = "lpSolveAPI")
#'
#' ## Estimate the bounds
#' obsEqMin(sset = sSet$sset,
#'          lpobj = lpObject,
#'          lpsolver = "lpSolveAPI")
#'
#' @export
obsEqMin <- function(sset, orig.sset = NULL, orig.criterion = NULL,
                     criterion.tol = 0, lpobj, lpsolver, lpsolver.options,
                     debug = FALSE) {
    if (!is.null(orig.sset)) {
        ## Prepare to obtain 'recentered' bootstrap
        ## criterion. Specifically, the |S| equality constraints are
        ## centered. Then, the original |S| equality constraints are
        ## added. In addition, 2 * |S| residual variables are added to
        ## the problem. These new residual variables correspond to the
        ## |S| equality constraints from the original, uncentered
        ## sample.
        tmpA <- NULL
        tmpRhs <- NULL
        tmpSense <- NULL
        scount <- 0
        for (s in names(orig.sset)) {
            avec <- replicate(2 * 2 * length(orig.sset), 0)
            avec[(2 * scount + 1):(2 * scount + 2)] <- c(-1, 1)
            avec <- c(avec, orig.sset[[s]]$g0, orig.sset[[s]]$g1)
            tmpA <- rbind(tmpA, avec)
            tmpRhs <- c(tmpRhs, orig.sset[[s]]$beta)
            tmpSense <- c(tmpSense, "=")
            scount <- scount + 1
        }
        avec <- c(rep(1, 2 * length(orig.sset)),
                  rep(0, 2 * length(orig.sset)),
                  rep(0, length(sset$s1$g0) + length(sset$s1$g1)))
        ## Update lpobj
        lpobj$ub <- c(rep(Inf, 2 * length(orig.sset)), lpobj$ub)
        lpobj$lb <- c(rep(0, 2 * length(orig.sset)), lpobj$lb)
        lpobj$A <- rbind(avec, tmpA,
                         cbind(matrix(0,
                                   nrow = nrow(lpobj$A),
                                   ncol = length(orig.sset) * 2),
                               lpobj$A))
        lpobj$rhs <- c(orig.criterion * (1 + criterion.tol),
                       tmpRhs, lpobj$rhs)
        lpobj$sense <- c("<=", tmpSense, lpobj$sense)
        lpobj$obj <- c(rep(0, 2 * length(orig.sset)), lpobj$obj)
    }
    lpsolver <- tolower(lpsolver)
    if (lpsolver == "gurobi") {
        if (debug & lpsolver.options$outputflag == 1) {
            message("\nMinimum criterion optimization statistics:")
            message("------------------------------------------")
        }
        model <- list()
        model$modelsense <- "min"
        model$obj   <- lpobj$obj
        model$A     <- lpobj$A
        lpobj$A <- NULL
        model$rhs   <- lpobj$rhs
        model$sense <- lpobj$sense
        model$ub    <- lpobj$ub
        model$lb    <- lpobj$lb
        if (debug) {
            gurobi::gurobi_write(model, "lpCriterion.mps")
            save(model, file = "lpCriterion.Rdata")
        }
        result   <- gurobi::gurobi(model, lpsolver.options)
        obseqmin <- result$objval
        optx     <- result$x
        status   <- result$status
        if (debug) cat("\n")
    }
    if (lpsolver == "cplexapi") {
        result <- runCplexAPI(lpobj, cplexAPI::CPX_MIN, lpsolver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
    }
    if (lpsolver == "lpsolveapi") {
        result <- runLpSolveAPI(lpobj, 'min', lpsolver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
    }
    ## provide nicer output
    g0sol <- optx[(2 * lpobj$sn + 1) : (2 * lpobj$sn + lpobj$gn0)]
    g1sol <- optx[(2 * lpobj$sn + lpobj$gn0 + 1) :
                  (2 * lpobj$sn + lpobj$gn0 + lpobj$gn1)]
    names(g0sol) <- names(sset$gstar$g0)
    names(g1sol) <- names(sset$gstar$g1)
    return(list(obj = obseqmin,
                g0 = g0sol,
                g1 = g1sol,
                status = status))
    ## object 'result' will not be returned---unnecessary, and very
    ## memory intensive.
}

#' Obtaining TE bounds
#'
#' This function estimates the bounds on the target treatment effect.
#' @param g0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param g1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set. This
#'     object is only used to determine the names of terms. If it is
#'     no submitted, then no names are provided to the solution
#'     vector.
#' @param lpobj A list of matrices and vectors defining an LP problem.
#' @param obseq.factor overall multiplicative factor for how much more
#'     the solution is permitted to violate observational equivalence
#'     of the IV-like estimands, i.e. \code{obseq.factor} will
#'     multiply \code{minobseq} directly.
#' @param noisy boolean, set to \code{TRUE} if optimization results
#'     should be displayed.
#' @param lpsolver string, name of the package used to solve the LP
#'     problem.
#' @param lpsolver.options list, each item of the list should
#'     correspond to an option specific to the LP solver
#'     selected.
#' @param debug boolean, indicates whether or not the function should
#'     provide output when obtaining bounds. The option is only
#'     applied when \code{lpsolver = 'gurobi'}. The output provided is
#'     the same as what the Gurobi API would send to the console.
#' @return a list containing the bounds on the treatment effect; the
#'     coefficients on each term in the MTR associated with the upper
#'     and lower bounds, for both counterfactuals; the optimization
#'     status to the maximization and minimization problems; the LP
#'     problem that the optimizer solved.
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#'
#' ## Declare empty list to be updated (in the event multiple IV like
#' ## specifications are provided
#' sSet <- list()
#'
#' ## Declare MTR formulas
#' formula1 = ~ 1 + u
#' formula0 = ~ 1 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#' polynomials1 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' ## Generate propensity score model
#' propensityObj <- propensity(formula = d ~ z,
#'                             data = dtm,
#'                             link = "linear")
#'
#' ## Generate IV estimates
#' ivEstimates <- ivEstimate(formula = ey ~ d | z,
#'                           data = dtm,
#'                           components = l(intercept, d),
#'                           treat = d,
#'                           list = FALSE)
#'
#' ## Generate target gamma moments
#' targetGamma <- genTarget(treat = "d",
#'                          m0 = ~ 1 + u,
#'                          m1 = ~ 1 + u,
#'                          target = "atu",
#'                          data = dtm,
#'                          splinesobj = splinesList,
#'                          pmodobj = propensityObj,
#'                          pm0 = polynomials0,
#'                          pm1 = polynomials1,
#'                          point = FALSE)
#'
#' ## Construct S-set. which contains the coefficients and weights
#' ## corresponding to various IV-like estimands
#' sSet <- genSSet(data = dtm,
#'                 sset = sSet,
#'                 sest = ivEstimates,
#'                 splinesobj = splinesList,
#'                 pmodobj = propensityObj$phat,
#'                 pm0 = polynomials0,
#'                 pm1 = polynomials1,
#'                 ncomponents = 2,
#'                 scount = 1,
#'                 yvar = "ey",
#'                 dvar = "d",
#'                 means = TRUE)
#'
#' ## Define additional upper- and lower-bound constraints for the LP
#' ## problem
#' A <- matrix(0, nrow = 22, ncol = 4)
#' A <- cbind(A, rbind(cbind(1, seq(0, 1, 0.1)),
#'                     matrix(0, nrow = 11, ncol = 2)))
#' A <- cbind(A, rbind(matrix(0, nrow = 11, ncol = 2),
#'                     cbind(1, seq(0, 1, 0.1))))
#'
#' sense <- c(rep(">", 11), rep("<", 11))
#' rhs <- c(rep(0.2, 11), rep(0.8, 11))
#'
#' ## Construct LP object to be interpreted and solved by lpSolveAPI
#' lpObject <- lpSetup(sset = sSet$sset,
#'                     mbA = A,
#'                     mbs = sense,
#'                     mbrhs = rhs,
#'                     lpsolver = "lpSolveAPI")
#'
#' ## Estimate the bounds
#' bound(g0 = targetGamma$gstar0,
#'       g1 = targetGamma$gstar1,
#'       sset = sSet$sset,
#'       lpobj = lpObject,
#'       obseq.factor = 1,
#'       lpsolver = "lpSolveAPI")
#'
#' @export
bound <- function(g0, g1, sset, lpobj, obseq.factor, lpsolver,
                  lpsolver.options, noisy = FALSE,
                  debug = FALSE) {
    lpsolver <- tolower(lpsolver)
    ## define model
    model <- list()
    model$obj <- c(replicate(2 * lpobj$sn, 0), g0, g1)
    model$rhs <- c(obseq.factor, lpobj$rhs)
    model$ub  <- lpobj$ub
    model$lb  <- lpobj$lb
    avec <- c(replicate(2 * lpobj$sn, 1),
              replicate(lpobj$gn0 + lpobj$gn1, 0))
    model$A <- rbind(avec, lpobj$A)
    model$sense <- c("<=", lpobj$sense)
    ## check scaling of model
    tmpA <- c(matrix(model$A, ncol = 1))
    tmpA <- tmpA[tmpA != 0]
    magDiffA <- max(magnitude(tmpA), na.rm = TRUE) -
        min(magnitude(tmpA), na.rm = TRUE)
    magDiffObj <- max(magnitude(model$obj), na.rm = TRUE) -
        min(magnitude(model$obj), na.rm = TRUE)
    magDiffRhs <- max(magnitude(model$rhs), na.rm = TRUE) -
        min(magnitude(model$rhs), na.rm = TRUE)
    modelStats  <-
        data.frame(matrix(c(min(abs(tmpA)),
                            min(abs(model$rhs[model$rhs != 0])),
                            min(abs(model$obj[model$obj != 0])),
                            max(abs(tmpA)),
                            max(abs(model$rhs[model$rhs != 0])),
                            max(abs(model$obj[model$obj != 0])),
                            magDiffA, magDiffRhs, magDiffObj),
                          ncol = 3))
    colnames(modelStats) <- c("Min abs.", "Max abs.", "Magnitude diff.")
    rownames(modelStats) <- c("Constraint matrix",
                              "RHS",
                              "Objective")
    rm(tmpA)
    ## obtain lower and upper bounds
    if (lpsolver == "gurobi") {
        if (debug & lpsolver.options$outputflag == 1) {
            message("\nLower bound optimization statistics:")
            message("------------------------------------")
        }
        if (debug == TRUE){
            model$modelsense <- "max"
            gurobi::gurobi_write(model, "lpMax.mps")
            save(model, file = "lpMax.Rdata")
            model$modelsense <- "min"
            gurobi::gurobi_write(model, "lpMin.mps")
            save(model, file = "lpMin.Rdata")
        }
        model$modelsense <- "min"
        minresult <- gurobi::gurobi(model, lpsolver.options)
        min <- minresult$objval
        minstatus <- 0
        if (minresult$status == "OPTIMAL") minstatus <- 1
        minoptx <- minresult$x
        if (debug & lpsolver.options$outputflag == 1) {
            message("\nUpper bound optimization statistics:")
            message("------------------------------------")
        }
        model$modelsense <- "max"
        maxresult <- gurobi::gurobi(model, lpsolver.options)
        max <- maxresult$objval
        maxstatus <- 0
        if (maxresult$status == "OPTIMAL") maxstatus <- 1
        maxoptx <- maxresult$x
        if (debug) cat("\n")
    }
    if (lpsolver == "cplexapi") {
        minresult <- runCplexAPI(model, cplexAPI::CPX_MIN, lpsolver.options)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runCplexAPI(model, cplexAPI::CPX_MAX, lpsolver.options)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
    }
    if (lpsolver == "lpsolveapi") {
        minresult <- runLpSolveAPI(model, 'min', lpsolver.options)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runLpSolveAPI(model, 'max', lpsolver.options)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
    }
    if (maxstatus == 0 || minstatus == 0) {
        return(NULL)
    }
    ming0 <- minoptx[(2 * lpobj$sn + 1) : (2 * lpobj$sn + lpobj$gn0)]
    ming1 <- minoptx[(2 * lpobj$sn + lpobj$gn0 + 1) :
                     (2 * lpobj$sn + lpobj$gn0 + lpobj$gn1)]
    maxg0 <- maxoptx[(2 * lpobj$sn + 1) : (2 * lpobj$sn + lpobj$gn0)]
    maxg1 <- maxoptx[(2 * lpobj$sn + lpobj$gn0 + 1) :
                     (2 * lpobj$sn + lpobj$gn0 + lpobj$gn1)]
    if (hasArg(sset)) {
        names(ming0) <- names(sset$s1$g0)
        names(ming1) <- names(sset$s1$g1)
        names(maxg0) <- names(sset$s1$g0)
        names(maxg1) <- names(sset$s1$g1)
    }
    if (noisy) {
        cat("Min status: ", minstatus, "\n", sep = "")
        cat("Max status: ", maxstatus, "\n", sep = "")
        cat("Bound: (", min, ", ", max, ")\n", sep = "")
    }
    return(list(max = max,
                maxg0 = maxg0,
                maxg1 = maxg1,
                maxresult = maxresult,
                maxstatus = maxstatus,
                min = min,
                ming0 = ming0,
                ming1 = ming1,
                minresult = minresult,
                minstatus = minstatus,
                model = model,
                modelstats = modelStats))
}


#' Running Gurobi LP solver
#'
#' This function solves the LP problem using the Gurobi package. The
#' object generated by \code{\link{lpSetup}} is compatible with the
#' \code{gurobi} function.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param lpsolver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from Gurobi. This includes the
#'     optimization status, the objective value, the solution vector,
#'     amongst other things.
runGurobi <- function(lpobj, lpsolver.options) {
    result <- gurobi::gurobi(lpobj, lpsolver.options)
    status <- 0
    if (result$status == "OPTIMAL") status <- 1
    optx <- result$x
    return(list(objval = result$objval,
                optx = result$x,
                status = status))
}


#' Running cplexAPI LP solver
#'
#' This function solves the LP problem using the cplexAPI package. The
#' object generated by \code{\link{lpSetup}} is not compatible
#' with the \code{cplexAPI} functions. This function adapts the object
#' to solve the LP problem.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param lpdir input either CPX_MAX or CPX_MIN, which sets the LP
#'     problem as a maximization or minimization problem.
#' @param lpsolver.options list, each item of the list should
#'     correspond to an option specific to the LP solver
#'     selected.
#' @return a list of the output from CPLEX. This includes the
#'     optimization status, the objective value, the solution vector,
#'     amongst other things.
runCplexAPI <- function(lpobj, lpdir, lpsolver.options) {
    ## Declare environment and set options
    env  <- cplexAPI::openEnvCPLEX()
    prob <- cplexAPI::initProbCPLEX(env)
    cplexAPI::chgProbNameCPLEX(env, prob, "sample")
    if (!is.null(lpsolver.options)) {
        for(i in seq(length(lpsolver.options))) {
            eval(parse(text = lpsolver.options[[i]]))
        }
    }
    ## Declare LP prblem
    sense <- lpobj$sense
    ## Original ------------------------
    ## Should no longer be needed since we now have the lpSetupSolver
    ## sense[sense == "<"]  <- "L"
    ## sense[sense == "<="] <- "L"
    ## sense[sense == ">"]  <- "G"
    ## sense[sense == ">="] <- "G"
    ## sense[sense == "="]  <- "E"
    ## sense[sense == "=="] <- "E"
    ## ub <- lpobj$ub
    ## ub[ub == Inf] <- cplexAPI::CPX_INFBOUND
    ## lb <- lpobj$lb
    ## lb[lb == -Inf] <- - cplexAPI::CPX_INFBOUND
    ## End original ---------------------------
    cnt <- apply(lpobj$A, MARGIN = 2, function(x) length(which(x != 0)))
    beg <- rep(0, ncol(lpobj$A))
    beg[-1] <- cumsum(cnt[-length(cnt)])
    ind <- unlist(apply(lpobj$A, MARGIN = 2, function(x) which(x != 0) - 1))
    print("Super inconvneinent input approahc for cplexapi")
    val <- c(lpobj$A)
    val <- val[val != 0]
    cplexAPI::copyLpwNamesCPLEX(env = env,
                                lp = prob,
                                nCols = ncol(lpobj$A),
                                nRows = nrow(lpobj$A),
                                lpdir = lpdir,
                                objf = lpobj$obj,
                                rhs = lpobj$rhs,
                                sense = sense,
                                matbeg = beg,
                                matcnt = cnt,
                                matind = ind,
                                matval = val,
                                lb = lb,
                                ub = ub)
    cplexAPI::lpoptCPLEX(env, prob)
    solution <- cplexAPI::solutionCPLEX(env, prob)
    cplexAPI::delProbCPLEX(env, prob)
    cplexAPI::closeEnvCPLEX(env)
    if (typeof(solution) == "S4") {
        if (attr(solution, "class") == "cplexError") {
            status <- 0
            solution <- list()
            solution$objval <- NA
            solution$x <- NA
        }
    }  else {
        if (solution$lpstat == 1) status <- 1
        if (solution$lpstat != 1) status <- 0
    }
    return(list(objval = solution$objval,
                optx   = solution$x,
                status = status))
}

#' Running lpSolveAPI
#'
#' This function solves the LP problem using the \code{lpSolveAPI}
#' package. The object generated by \code{\link{lpSetup}} is not
#' compatible with the \code{lpSolveAPI} functions. This function
#' adapts the object to solve the LP problem.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param modelsense input either 'max' or 'min' which sets the LP
#'     problem as a maximization or minimization problem.
#' @param lpsolver.options list, each item of the list should
#'     correspond to an option specific to the LP solver
#'     selected.
#' @return a list of the output from \code{lpSolveAPI}. This includes
#'     the optimization status, the objective value, the solution
#'     vector.
runLpSolveAPI <- function(lpobj, modelsense, lpsolver.options) {
    lpmodel <- lpSolveAPI::make.lp(nrow(lpobj$A), ncol(lpobj$A))
    for (j in 1:ncol(lpobj$A)) {
        lpSolveAPI::set.column(lprec = lpmodel,
                               column = j,
                               x = lpobj$A[, j])
    }
    lpSolveAPI::set.constr.value(lprec = lpmodel,
                                 rhs = lpobj$rhs)
    sense <- lpobj$sense
    sense[sense == "<"]  <- "<="
    sense[sense == ">"]  <- ">="
    sense[sense == "=="] <- "="
    lpSolveAPI::set.constr.type(lprec = lpmodel,
                                types = sense)
    lpSolveAPI::set.objfn(lprec = lpmodel,
                          obj = lpobj$obj)
    lpSolveAPI::lp.control(lprec = lpmodel,
                           sense = modelsense)
    if (!is.null(lpsolver.options)) {
        eval(lpsolver.options)
    }
    lpSolveAPI::set.bounds(lprec = lpmodel,
                           lower = lpobj$lb,
                           upper = lpobj$ub)
    solved <- lpSolveAPI::solve.lpExtPtr(lpmodel)
    if (solved == 0) status <- 1
    if (solved != 0) status <- 0
    return(list(objval = lpSolveAPI::get.objective(lpmodel),
                optx   = lpSolveAPI::get.variables(lpmodel),
                status = status))
}

#' Check magnitude of real number
#'
#' This function returns the order of magnitude of a a number.
#'
#' @param x The number to be checked.
#' @return An integer indicating the order of magnitude.
magnitude <- function(x) {
    sapply(x, function(y) {
        if (y == 0) return(NA)
        else return(floor(log(abs(y), 10)))
    })
}

#' Function to parse options for Gurobi
#'
#' This function constructs a list of options to be parsed when
#' \code{lpsolver} is set to \code{Gurobi}. This function really
#' implements some default values, and accounts for the \code{debug}
#' option.
#' @param options list. The list should be structured the same way as
#'     if one were using the \code{gurobi} library directly. That is,
#'     the name of each item must be the name of the option, and is
#'     case sensitive. The value assigned to each item is the value to
#'     set the option to.
#' @return list, the set of options declared by the user, including
#'     some additional default values (if not assigned by the user)
#'     and accounting for \code{debug}.
optionsGurobi <- function(options, debug) {
    if (! "outputflag" %in% names(options)) {
        if (debug)  options$outputflag = 1
        if (!debug) options$outputflag = 0
    }
    if (! "dualreductions" %in% names(options)) {
        options$dualreductions <- 1
    }
    if (! "FeasibilityTol" %in% names(options)) {
        options$FeasibilityTol <- 1e-06
    }
    if (! "presolve" %in% names(options)) {
        if (hasArg(lpsolver.presolve)) {
            options$presolve <- as.integer(lpsolver.presolve)
        }
    }
    return(options)
}

#' Function to parse options for lp_solve
#'
#' This function constructs a list of options to be parsed when
#' \code{lpsolver} is set to \code{lpsolveapi}. The options permitted
#' are those that can be set via \code{lpSolveAPI::lp.control}, and
#' should be passed as a named list (e.g. \code{list(epslevel =
#' "tight")}).
#' @param options list. The name of each item must be the name of the
#'     option, and is case sensitive. The value assigned to each item
#'     is the value to set the option to. The \code{lprec} argument
#'     should always be omitted.
#' @return string, the command to be evaluated to implement the
#'     options.
optionsLpSolveAPI <- function(options) {
    ## Implement default tolerance
    if (!"epslevel" %in% names(options)) {
        options$epslevel = "tight"
    }
    optionsStr <- gsub("\\s+", " ", Reduce(paste, deparse(options)))
    optionsStr <- gsub("list\\(", "lpSolveAPI::lp.control(lprec = lpmodel, ",
                       optionsStr)
    return(optionsStr)
}

#' Function to parse options for CPLEX
#'
#' This function constructs a list of options to be parsed when
#' \code{lpsolver} is set to \code{cplexapi}.
#' @param options list. The name of each item must be the name of the
#'     function to set the option, and is case sensitive. The value
#'     assigned to each item is the value to set the option to. The
#'     \code{env} argument should always be omitted. If the option
#'     accepts a list of parameters, then these parameters should be
#'     passed as using a named vector (e.g.
#'     \code{list(setLogFileNameCPLEX = c(filename = "cpx.log", mode =
#'     "w"))}).  If the function to set the option can be used
#'     multiple times, then the value submitted should be a a list,
#'     with each entry being a named vector
#'     (e.g. \code{list(setDblParmCPLEX = list(c(parm = 1016, value =
#'     1e-04), c(parm = 1084, value = 2)))}). If the option only
#'     requires the \code{env} parameter, then an \code{NA} should be
#'     passsed as the parameter value (e.g. \code{list(setDefaultParm
#'     = NA)}).
#' @return list, each element being the command to evaluate to
#'     implement an option.
optionsCplexAPI <- function(options) {
    ## Implement default tolerance
    if ("setDblParmCPLEX" %in% names(options)) {
        pos <- which(names(options) == "setDblParmCPLEX")
        if (!is.list(options[[pos]])) {
            if (options[[pos]]["parm"] != 1016) {
                options[[pos]] <- list(options[[pos]],
                                       c(parm = 1016, value = 1e-06))
            }
        } else {
            parms <- unlist(lapply(options[[pos]], function(x) x["parm"]))
            if (! 1016 %in% parms) {
                options[[pos]][length(options[[pos]]) + 1] <-
                    c(parm = 1016, value = 1e-06)
            }
        }
    } else {
        options$setDblParmCPLEX <- c(parm = 1016, value = 1e-06)
    }
    ## Construct commands to implement options
    optionsStr <- list()
    counter <- 1
    for (i in seq(length(options))) {
        if (!is.list(options[[i]])) {
            optionsStr[counter] <-
                optionsCplexAPISingle(names(options)[i],
                                      options[[i]])
            counter <- counter + 1
        }
        if (is.list(options[[i]])) {
            for (j in seq(length(options[[i]]))) {
                optionsStr[counter] <-
                    optionsCplexAPISingle(names(options)[i],
                                          options[[i]][[j]])
                counter <- counter + 1
            }
        }
    }
    return(optionsStr)
}

#' Function to parse a single set of options for CPLEX
#'
#' This function constructs a string to be parsed when \code{lpsolver}
#' is set to \code{cplexapi}.
#' @param name string, name of the \code{cplexapi} function to call to
#'     implement the option.
#' @param vector a named vector, contains the argument names and
#'     values of the options. The \code{env} argument in the
#'     \code{cplexapi} documentation should always be omitted.
#' @return string, the command to be evaluated to implement a single
#'     option.
optionsCplexAPISingle <- function(name, vector) {
    suppressWarnings(
    if (is.na(vector)) {
        tmpCommand <- paste0("cplexAPI::",
                             name,
                             "(env = env)")
    } else {
        tmpCommand <- NULL
        for (i in seq(length(vector))) {
            if (!is.character(vector[i])) {
                tmpAdd <- paste0(names(vector)[i], " = ", vector[[i]])
            } else {
                tmpAdd <- paste0(names(vector)[i], " = ", deparse(vector[[i]]))
            }
            if (i == 1) tmpCommand <- tmpAdd
            if (i != 1) tmpCommand <- paste(tmpCommand, tmpAdd, sep = ", ")
        }
        tmpCommand <- paste0("cplexAPI::",
                             name,
                             "(env = env, ",
                             tmpCommand,
                             ")")
    }
    )
    return(tmpCommand)
}

#' Function to extract feasibility tolerance from CPLEX options
#'
#' This function parses through the user-submitted CPLEX options to
#' determine what the feasibility tolerance is. This tolerance can
#' then be used for the audit.  If the user does not set the CPLEX
#' feasibility tolerance, then a default value of \code{1e-06} is
#' returned.
#' @param options list, the set of options submitted by the user.
#' @return scalar, the level to set the audit tolerance at.
optionsCplexAPITol <- function(options) {
    if (! "setDblParmCPLEX" %in% names(options)) {
        audit.tol <- 1e-06
    } else {
        if (is.list(options$setDblParmCPLEX)) {
            parms <- unlist(lapply(options$setDblParmCPLEX,
                                   function(x) x["parm"]))
            if (1016 %in% parms) {
                pos <- which(parms == 1016)
                audit.tol <- options$setDblParmCPLEX[[pos]]["value"]
            } else {
                audit.tol <- 1e-06
            }
        } else {
            if (options$setDblParmCPLEX["parm"] == 1016) {
                audit.tol <- options$setDblParmCPLEX["value"]
            } else {
                audit.tol <- 1e-06
            }
        }
    }
    return(audit.tol)
}
