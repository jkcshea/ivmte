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
        A     <- rbind(A, mbA)
        sense <- c(sense, mbs)
        rhs   <- c(rhs, mbrhs)
    }
    ## Adjust for Rcplex
    if (lpsolver == "rcplex") {
        sense[sense == "<"]  <- "L"
        sense[sense == "<="] <- "L"
        sense[sense == ">"]  <- "G"
        sense[sense == ">="] <- "G"
        sense[sense == "="]  <- "E"
    }
    if (lpsolver %in% c("gurobi", "rcplex", "lpsolveapi")) {
        A <- Matrix::Matrix(A, sparse = TRUE)
    }
    return(list(obj = obj,
                rhs = rhs,
                sense = sense,
                A = A,
                ub = ub,
                lb = lb,
                sn = sn,
                gn0 = gn0,
                gn1 = gn1))
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
#' @param obseq.tol tolerance for violation of observational
#'     equivalence, set to 0 by default.
#' @param lpobj A list of matrices and vectors defining an LP problem.
#' @param lpsolver string, name of the package used to solve the LP
#'     problem.
#' @return A list including the minimum violation of observational
#'     equivalence, the solution to the LP problem, and the status of
#'     the solution.
#'
#' @examples
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
#'                          uname = u,
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
                     obseq.tol = 0, lpobj, lpsolver) {
    if (!is.null(orig.sset)) {
        ## Prepare to obtain 'recentered' bootstrap criterion
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
        lpobj$rhs <- c(orig.criterion * (1 + obseq.tol),
                       tmpRhs, lpobj$rhs)
        lpobj$sense <- c("<=", tmpSense, lpobj$sense)
        lpobj$obj <- c(rep(0, 2 * length(orig.sset)), lpobj$obj)
    }
    lpsolver <- tolower(lpsolver)
    if (lpsolver == "gurobi") {
        model <- list()
        model$modelsense <- "min"
        model$obj   <- lpobj$obj
        model$A     <- lpobj$A
        model$rhs   <- lpobj$rhs
        model$sense <- lpobj$sense
        model$ub    <- lpobj$ub
        model$lb    <- lpobj$lb
        result   <- gurobi::gurobi(model, list(outputflag = 0))
        obseqmin <- result$objval
        optx     <- result$x
        status   <- result$status
    }
    if (lpsolver == "rcplex") {
        result <- Rcplex::Rcplex(objsense = "min",
                                 cvec = lpobj$obj,
                                 Amat = lpobj$A,
                                 bvec = lpobj$rhs,
                                 sense = lpobj$sense,
                                 ub = lpobj$ub,
                                 lb = lpobj$lb,
                                 control = list(trace = FALSE))
        obseqmin <- result$obj
        optx     <- result$xopt
        status   <- result$status
    }
    if (lpsolver == "cplexapi") {
        result <- runCplexAPI(lpobj, cplexAPI::CPX_MIN)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
    }
    if (lpsolver == "lpsolveapi") {
        result <- runLpSolveAPI(lpobj, 'min')
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
                status = status,
                result = result))
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
#' @return a list containing the bounds on the treatment effect; the
#'     coefficients on each term in the MTR associated with the upper
#'     and lower bounds, for both counterfactuals; the optimization
#'     status to the maximization and minimization problems; the LP
#'     problem that the optimizer solved.
#'
#' @examples
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
#'                          uname = u,
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
bound <- function(g0, g1, sset, lpobj, obseq.factor, lpsolver, noisy = FALSE) {
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
    ## obtain lower and upper bounds
    if (lpsolver == "gurobi") {
        model$modelsense <- "min"
        minresult <- gurobi::gurobi(model, list(outputflag = 0,
                                                dualreductions = 1))
        min <- minresult$objval
        minstatus <- 0
        if (minresult$status == "OPTIMAL") minstatus <- 1
        minoptx <- minresult$x
        model$modelsense <- "max"
        maxresult <- gurobi::gurobi(model, list(outputflag = 0,
                                                dualreductions = 1))
        max <- maxresult$objval
        maxstatus <- 0
        if (maxresult$status == "OPTIMAL") maxstatus <- 1
        maxoptx <- maxresult$x
    }
    if (lpsolver == "rcplex") {
        model$sense[1] <- "L" ## to satisfy minimal obs. equiv. deviation
        minresult <- Rcplex::Rcplex(objsense = "min",
                                    cvec = model$obj,
                                    Amat = model$A,
                                    bvec = model$rhs,
                                    sense = model$sense,
                                    ub = model$ub,
                                    lb = model$lb,
                                    control = list(trace = FALSE))
        min <- minresult$obj
        minstatus <- minresult$status
        minoptx   <- minresult$xopt
        maxresult <- Rcplex::Rcplex(objsense = "max",
                                    cvec = model$obj,
                                    Amat = model$A,
                                    bvec = model$rhs,
                                    sense = model$sense,
                                    ub = model$ub,
                                    lb = model$lb,
                                    control = list(trace = FALSE))
        max <- maxresult$obj
        maxstatus <- maxresult$status
        maxoptx   <- maxresult$xopt
        ## Deal with cases where the bounds are contradictory,
        ## which can occur when the target parameter is
        ## unbounded. However, this can potentially be a precision
        ## issue, i.e. the min can exceed the max by a trivial
        ## amount. You suspect this is a threshold problem, so you
        ## allow for some tolerance.
        if (min - max > 0) {
            if (round(min - max, 10) > 0) {
                min <- NULL
                max <- NULL
                minstatus <- 0
                maxstatus <- 0
                minoptx <- NULL
                maxoptx <- NULL
            } else {
                max <- min
            }
        }
    }
    if (lpsolver == "cplexapi") {
        minresult <- runCplexAPI(model, cplexAPI::CPX_MIN)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runCplexAPI(model, cplexAPI::CPX_MAX)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
    }
    if (lpsolver == "lpsolveapi") {
        minresult <- runLpSolveAPI(model, 'min')
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runLpSolveAPI(model, 'max')
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
                model = model))
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
#' @return a list of the output from CPLEX. This includes the
#'     optimization status, the objective value, the solution vector,
#'     amongst other things.
runCplexAPI <- function(lpobj, lpdir) {
    env  <- cplexAPI::openEnvCPLEX()
    prob <- cplexAPI::initProbCPLEX(env)
    cplexAPI::chgProbNameCPLEX(env, prob, "sample")
    sense <- lpobj$sense
    sense[sense == "<"]  <- "L"
    sense[sense == "<="] <- "L"
    sense[sense == ">"]  <- "G"
    sense[sense == ">="] <- "G"
    sense[sense == "="]  <- "E"
    sense[sense == "=="] <- "E"
    ub <- lpobj$ub
    ub[ub == Inf] <- cplexAPI::CPX_INFBOUND
    lb <- lpobj$lb
    lb[lb == -Inf] <- - cplexAPI::CPX_INFBOUND
    cnt <- apply(lpobj$A, MARGIN = 2, function(x) length(which(x != 0)))
    beg <- rep(0, ncol(lpobj$A))
    beg[-1] <- cumsum(cnt[-length(cnt)])
    ind <- unlist(apply(lpobj$A, MARGIN = 2, function(x) which(x != 0) - 1))
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
#' @return a list of the output from \code{lpSolveAPI}. This includes
#'     the optimization status, the objective value, the solution
#'     vector.
runLpSolveAPI <- function(lpobj, modelsense) {
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
