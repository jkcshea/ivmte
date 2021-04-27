#' Constructing LP problem
#'
#' This function takes in the IV estimates from the set of IV
#' regressions declared by the user, as well as their corresponding
#' moments of the terms in the MTR. These are then used to construct
#' the components that make up the LP problem. Note that the LP model
#' will be saved inside an environment variable, which is to be passed
#' through the argument \code{env}. This is done for efficient use of
#' memory. The environment \code{env} is supposed to already contain a
#' list under the entry \code{$mbobj} containing the matrices defining
#' the shape constraints. This list of shape constraints \code{$mbobj}
#' should contain three entries corresponding to a system of linear
#' equations of the form \code{Ax <=> b}: \code{mbA}, the matrix
#' defining the constraints, \code{A}; \code{mbs}, a vector indicating
#' whether a row in \code{mbA} is an equality or inequality constraint
#' (for Gurobi and CPLEX, use '<=', '>=', '='; for lpSolveAPI, use
#' 'L', 'G', and 'E'); \code{mbrhs}, a vector of the right hand side
#' values defining the constraint of the form i.e. the vector
#' \code{b}. Depending on the linear programming solver used, this
#' function will return different output specific to the solver.
#' @param env environment containing the matrices defining the LP
#'     problem.
#' @param sset List of IV-like estimates and the corresponding gamma
#'     terms.
#' @param orig.sset list, only used for bootstraps. The list contains
#'     the gamma moments for each element in the S-set, as well as the
#'     IV-like coefficients.
#' @param shape boolean, default set to TRUE. Switch to determine
#'     whether or not to include shape restrictions in the LP problem.
#' @param direct boolean, set to \code{TRUE} if the direct MTR
#'     regression is used.
#' @param rescale boolean, set to \code{TRUE} if the MTR components
#'     should be rescaled to improve stability in the LP/QP/QCP
#'     optimization.
#' @param solver string, name of the package used to solve the LP
#'     problem.
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
#' formula0 = ~ 1 + u
#' formula1 = ~ 1 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                           data = dtm,
#'                           uname = u,
#'                           as.function = FALSE)
#' polynomials1 <- polyparse(formula = formula1,
#'                           data = dtm,
#'                           uname = u,
#'                            as.function = FALSE)
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
#' ## Only the entry $sset is required
#' sSet <- sSet$sset
#'
#' ## Define additional upper- and lower-bound constraints for the LP
#' ## problem.  The code below imposes a lower bound of 0.2 and upper
#' ## bound of 0.8 on the MTRs.
#' A <- matrix(0, nrow = 22, ncol = 4)
#' A <- cbind(A, rbind(cbind(1, seq(0, 1, 0.1)),
#'                     matrix(0, nrow = 11, ncol = 2)))
#' A <- cbind(A, rbind(matrix(0, nrow = 11, ncol = 2),
#'                     cbind(1, seq(0, 1, 0.1))))
#' sense <- c(rep(">", 11), rep("<", 11))
#' rhs <- c(rep(0.2, 11), rep(0.8, 11))
#'
#' ## Construct LP object to be interpreted and solved by
#' ## lpSolveAPI. Note that an environment has to be created for the LP
#' ## object. The matrices defining the shape restrictions must be stored
#' ## as a list under the entry \code{$mbobj} in the environment.
#' lpEnv <- new.env()
#' lpEnv$mbobj <- list(mbA = A,
#'                     mbs = sense,
#'                     mbrhs = rhs)
#' ## Convert the matrices defining the shape constraints into a format
#' ## that is suitable for the LP solver.
#' lpSetup(env = lpEnv,
#'         sset = sSet,
#'         solver = "lpsolveapi")
#' ## Setup LP model so that it is solving for the bounds.
#' lpSetupBound(env = lpEnv,
#'              g0 = targetGamma$gstar0,
#'              g1 = targetGamma$gstar1,
#'              sset = sSet,
#'              criterion.tol = 0,
#'              criterion.min = 0,
#'              solver = "lpsolveapi")
#' ## Declare any LP solver options as a list.
#' lpOptions <- optionsLpSolveAPI(list(epslevel = "tight"))
#' ## Obtain the bounds.
#' bounds <- bound(env = lpEnv,
#'                 sset = sSet,
#'                 solver = "lpsolveapi",
#'                 solver.options = lpOptions)
#' cat("The bounds are [",  bounds$min, ",", bounds$max, "].\n")
#'
#' @export
lpSetup <- function(env, sset, orig.sset = NULL,
                    equal.coef0 = NULL, equal.coef1 = NULL,
                    shape = TRUE, direct = FALSE, rescale = TRUE,
                    solver) {
    ## Read in constraint grids and sequences
    solver <- tolower(solver)
    for (i in names(env$shapeSeq)) {
        if (length(env$shapeSeq[[i]]) == 0) {
            env$shapeSeq[[i]] <- NULL
        }
    }
    if (!direct) {
        ## Determine lengths
        sn  <- length(sset)
        gn0 <- length(sset$s1$g0)
        gn1 <- length(sset$s1$g1)
        ## Generate all vectors/matrices for LP optimization to minimize
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
        ## Generate informative colnumn names and row names indicating
        ## which rows correspond to which IV-like specifications, and
        ## which columns correspond to which MTR terms
        ## colnames(A) <- c(c(rbind(paste0('slack', seq(sn), '-'),
        ##                          paste0('slack', seq(sn), '+'))),
        ##                  colnames(A)[(2 * sn + 1) : ncol(A)])
        tmpIvs <- paste0('iv', lapply(sset, function(x) x$ivspec))
        tmpBetas <- lapply(sset, function(x) names(x$beta))
        rownames(A) <- mapply(paste, tmpIvs, tmpBetas, sep = '.')
        rm(tmpIvs, tmpBetas)
    } else {
        sn <- 0
        A <- NULL
        sense <- NULL
        rhs <- NULL
        gn0 <- ncol(sset$s1$g0)
        gn1 <- ncol(sset$s1$g1)
    }
    ## Add additional equality constraints if included
    if (!is.null(equal.coef0) & !is.null(equal.coef1)) {
        equal.coef0 <- paste0('[m0]', equal.coef0)
        equal.coef1 <- paste0('[m1]', equal.coef1)
        if (!direct) {
            tmpANames <- colnames(A)
        } else {
            tmpANames <- c(colnames(sset$s1$g0), colnames(sset$s1$g1))
        }
        equal.constraints <- lpSetupEqualCoef(equal.coef0, equal.coef1,
                                              tmpANames)
        A <- rbind(A, equal.constraints$A)
        sense <- c(sense, equal.constraints$sense)
        rhs <- c(rhs, equal.constraints$rhs)
        rm(tmpANames)
    }
    ## Add in additional constraints if included
    if (shape == TRUE) {
        mbA <- rbind(A, env$mbobj$mbA)
        if (direct) colnames(mbA) <- colnames(env$mbobj$mbA)
        env$mbobj$mbA <- NULL
        sense <- c(sense, env$mbobj$mbs)
        env$mbobj$mbs <- NULL
        rhs <- c(rhs, env$mbobj$mbrhs)
        env$mbobj$mbrhs <- NULL
    } else {
        mbA <- A
    }
    rm(A)
    if (!direct) {
        colnames(mbA) <- c(c(rbind(paste0('slack', seq(sn), '-'),
                                   paste0('slack', seq(sn), '+'))),
                           names(sset$s1$g0), names(sset$s1$g1))
    } else {
        colnames(mbA) <- c(colnames(sset$s1$g0), colnames(sset$s1$g1))
    }
    ## Define bounds on parameters
    ub <- replicate(ncol(mbA), Inf)
    lb <- c(unlist(replicate(sn * 2, 0)), replicate(gn0 + gn1, -Inf))
    if (direct && rescale) {
        ## Rescale linear constraints
        colNorms0 <- apply(sset$s1$g0, MARGIN = 2, function(x) sqrt(sum(x^2)))
        colNorms1 <- apply(sset$s1$g1, MARGIN = 2, function(x) sqrt(sum(x^2)))
        colNorms <- c(colNorms0, colNorms1)
        rm(colNorms0, colNorms1)
        mbA <- sweep(x = mbA,
                     MARGIN = 2,
                     STATS = colNorms,
                     FUN = '/')
        env$colNorms <- colNorms
    }
    ## Convert into sparse matrix
    if (solver %in% c("gurobi", "lpsolveapi")) {
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

#' Generate equality constraints
#'
#' This function generates the linear constraints to ensure that
#' certain MTR coefficients are constant across the treatment and
#' control group.
#' @param equal.coef0 character, name of terms in \code{m0} that
#'     should have common coefficients with the corresponding terms in
#'     \code{m1}.
#' @param equal.coef1 character, name of terms in \code{m1} that
#'     should have common coefficients with the corresponding terms in
#'     \code{m0}.
#' @param ANames character, name of all terms in \code{m0} and
#'     \code{m1}. The names of the terms corresponding to the
#'     treatment and control groups should be distinguishable. For
#'     example, all terms for \code{m0} may contain a prefix '[m0]',
#'     and all terms for \code{m1} may contain a prefix '[m1]'. All
#'     the terms in \code{equal.coef0} and \code{equal.coef1} should
#'     be contained in \code{ANames}.
#' @return A list, containing the matrix of linear equality
#'     constraints, a vector of equal signs, and a vector of 0s.
lpSetupEqualCoef <- function(equal.coef0, equal.coef1, ANames) {
    tmpPos0 <- sapply(equal.coef0, function(x) which(ANames == x))
    tmpPos1 <- sapply(equal.coef1, function(x) which(ANames == x))
    newA <- NULL
    newSense <- NULL
    newRhs <- NULL
    for (i in 1:length(tmpPos0)) {
        tmpA <- matrix(0, ncol = length(ANames), nrow = 1)
        tmpA[tmpPos0[i]] <- 1
        tmpA[tmpPos1[i]] <- -1
        rownames(tmpA) <- paste0('eq.coef.', i)
        newA <- rbind(newA, tmpA)
        newSense <- c(newSense, '=')
        newRhs <- c(newRhs, 0)
    }
    return(list(A = newA,
                sense = newSense,
                rhs = newRhs))
}

#' Configure LP environment for diagnostics
#'
#' This function separates the shape constraints from the LP
#' environment. That way, the model can be solved without any shape
#' constraints, which is the primary cause of infeasibility. This is
#' done in order to check which shape constraints are causing the
#' model to be infeasible. The LP model must be passed as an
#' environment variable, under the entry \code{$lpobj}. See
#' \code{\link{lpSetup}}.
#' @param env The LP environment
#' @param sset List of IV-like estimates and the corresponding gamma
#'     terms.
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
lpSetupInfeasible <- function(env, sset) {
    sn <- length(sset)
    ## Separate shape constraint objects
    env$mbobj$mbA <- env$lpobj$A[-(1:sn), ]
    env$mbobj$mbrhs <- env$lpobj$rhs[-(1:sn)]
    env$mbobj$mbsense <- env$lpobj$rhs[-(1:sn)]
    ## Reduce lpobjects so they do not contain any shape constraints
    env$lpobj$A <- env$lpobj$A[1:sn, ]
    env$lpobj$rhs <- env$lpobj$rhs[1:sn]
    env$lpobj$sense <- env$lpobj$sense[1:sn]
}

#' Configure LP environment for minimizing the criterion
#'
#' This function sets up the objective function for minimizing the
#' criterion. The LP model must be passed as an environment variable,
#' under the entry \code{$lpobj}. See \code{\link{lpSetup}}.
#' @param env The LP environment
#' @param sset List of IV-like estimates and the corresponding gamma
#'     terms.
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
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

#' Configure LP environment to be compatible with solvers
#'
#' This alters the LP environment so the model will be compatible with
#' specific solvers. The LP model must be passed as an environment
#' variable, under the entry \code{$lpobj}. See \code{\link{lpSetup}}.
#' @param env The LP environment
#' @param solver Character, the LP solver.
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
lpSetupSolver <- function(env, solver) {
    if (solver == "cplexapi") {
        env$lpobj$sense[env$lpobj$sense == "<"]  <- "L"
        env$lpobj$sense[env$lpobj$sense == "<="] <- "L"
        env$lpobj$sense[env$lpobj$sense == ">"]  <- "G"
        env$lpobj$sense[env$lpobj$sense == ">="] <- "G"
        env$lpobj$sense[env$lpobj$sense == "="]  <- "E"
        env$lpobj$sense[env$lpobj$sense == "=="] <- "E"
        env$lpobj$ub[env$lpobj$ub == Inf] <- cplexAPI::CPX_INFBOUND
        env$lpobj$lb[env$lpobj$lb == -Inf] <- -cplexAPI::CPX_INFBOUND
    }
    if (solver == "lpsolveapi") {
        env$lpobj$sense[env$lpobj$sense == "<"]  <- "<="
        env$lpobj$sense[env$lpobj$sense == ">"]  <- ">="
        env$lpobj$sense[env$lpobj$sense == "=="] <- "="
    }
}

#' Configure LP environment for specification testing
#'
#' This function re-centers various objects in the LP environment so
#' that a specification test can be performed via the bootstrap. The
#' LP model must be passed as an environment variable, under the entry
#' \code{$lpobj}. See \code{\link{lpSetup}}.
#' @param env the LP environment
#' @param sset list of IV-like estimates and the corresponding gamma
#'     terms.
#' @param orig.sset list, only used for bootstraps. The list caontains
#'     the gamma moments for each element in the S-set, as well as the
#'     IV-like coefficients.
#' @param orig.criterion scalar, only used for bootstraps. This is the
#'     minimum criterion from the original sample.
#' @param criterion.tol tolerance for violation of observational
#'     equivalence, set to 0 by default.
#' @param setup boolean. If \code{TRUE}, the function will modify the
#'     LP environment so that the LP solver can obtain the test
#'     statistic for the specification test. If \code{FALSE}, then it
#'     will undo the changes made by the function if \code{setup =
#'     TRUE}.
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
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


#' Configure LP environment for obtaining the bounds
#'
#' This function sets up the LP model so that the bounds can be
#' obtained. The LP model must be passed as an environment variable,
#' under the entry \code{$lpobj}. See \code{\link{lpSetup}}.
#' @param env the environment containing the LP model.
#' @param g0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param g1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set. This
#'     object is only used to determine the names of terms. If it is
#'     no submitted, then no names are provided to the solution
#'     vector.
#' @param criterion.tol additional multiplicative factor for how much
#'     more the solution is permitted to violate observational
#'     equivalence of the IV-like estimands, i.e. \code{1 +
#'     criterion.tol} will multiply \code{criterion.min} directly.
#' @param criterion.min minimum criterion, i.e. minimum deviation from
#'     observational equivalence while satisfying shape constraints.
#' @param solver string, name of the package used to solve the LP
#'     problem.
#' @param setup boolean. If \code{TRUE}, the function will modify the
#'     LP environment so that the LP solver can obtain the bounds. If
#'     \code{FALSE}, then it will undo the changes made by the
#'     function if \code{setup = TRUE}.
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
lpSetupBound <- function(env, g0, g1, sset, criterion.tol, criterion.min,
                         solver, setup = TRUE) {
    if (setup) {
        solver <- tolower(solver)
        ## Update objective function
        tmpSlack <- replicate(2 * env$lpobj$sn, 0)
        names(tmpSlack) <- c(rbind(paste0('slack', seq(env$lpobj$sn), '-'),
                                   paste0('slack', seq(env$lpobj$sn), '+')))
        env$lpobj$obj <- c(tmpSlack, g0, g1)
        ## Allow for slack in minimum criterion
        env$lpobj$rhs <- c(criterion.min * (1 + criterion.tol), env$lpobj$rhs)
        avec <- c(replicate(2 * env$lpobj$sn, 1),
                  replicate(env$lpobj$gn0 + env$lpobj$gn1, 0))
        env$lpobj$A <- rbind(avec, env$lpobj$A)
        ## Label each row with the corresponding constraint. IV-like
        ## constraints are already labeled.
        rownames(env$lpobj$A)[1] <- 'criterion'
        tmpOffset <- env$lpobj$sn + 1
        if (length(env$mbobj$lb0seq) > 0) {
            rownames(env$lpobj$A)[env$mbobj$lb0seq + tmpOffset] <- 'm0.lb'
        }
        if (length(env$mbobj$lb1seq) > 0) {
            rownames(env$lpobj$A)[env$mbobj$lb1seq + tmpOffset] <- 'm1.lb'
        }
        if (length(env$mbobj$lbteseq) > 0) {
            rownames(env$lpobj$A)[env$mbobj$lbteseq + tmpOffset] <- 'mte.lb'
        }
        if (length(env$mbobj$ub0seq) > 0) {
            rownames(env$lpobj$A)[env$mbobj$ub0seq + tmpOffset] <- 'm0.ub'
        }
        if (length(env$mbobj$ub1seq) > 0) {
            rownames(env$lpobj$A)[env$mbobj$ub1seq + tmpOffset] <- 'm1.ub'
        }
        if (length(env$mbobj$ubteseq) > 0) {
            rownames(env$lpobj$A)[env$mbobj$ubteseq + tmpOffset] <- 'mte.ub'
        }
        if (!is.null(env$mbobj$mono0seq)) {
            tmpDec <- which(env$mbobj$mono0seq[, 2] == -1)
            if (length(tmpDec) > 0) {
                tmpDec <- env$mbobj$mono0seq[, 1][tmpDec]
                rownames(env$lpobj$A)[tmpDec + tmpOffset] <- 'm0.dec'
            }
            rm(tmpDec)
            tmpInc <- which(env$mbobj$mono0seq[, 2] == 1)
            if (length(tmpInc) > 0) {
                tmpInc <- env$mbobj$mono0seq[, 1][tmpInc]
                rownames(env$lpobj$A)[tmpInc + tmpOffset] <- 'm0.inc'
            }
            rm(tmpInc)
        }
        if (!is.null(env$mbobj$mono1seq)) {
            tmpDec <- which(env$mbobj$mono1seq[, 2] == -1)
            if (length(tmpDec) > 0) {
                tmpDec <- env$mbobj$mono1seq[, 1][tmpDec]
                rownames(env$lpobj$A)[tmpDec + tmpOffset] <- 'm1.dec'
            }
            rm(tmpDec)
            tmpInc <- which(env$mbobj$mono1seq[, 2] == 1)
            if (length(tmpInc) > 0) {
                tmpInc <- env$mbobj$mono1seq[, 1][tmpInc]
                rownames(env$lpobj$A)[tmpInc + tmpOffset] <- 'm1.inc'
            }
            rm(tmpInc)
        }
        if (!is.null(env$mbobj$monoteseq)) {
            tmpDec <- which(env$mbobj$monoteseq[, 2] == -1)
            if (length(tmpDec) > 0) {
                tmpDec <- env$mbobj$monoteseq[, 1][tmpDec]
                rownames(env$lpobj$A)[tmpDec + tmpOffset] <- 'mte.dec'
            }
            rm(tmpDec)
            tmpInc <- which(env$mbobj$monoteseq[, 2] == 1)
            if (length(tmpInc) > 0) {
                tmpInc <- env$mbobj$monoteseq[, 1][tmpInc]
                rownames(env$lpobj$A)[tmpInc + tmpOffset] <- 'mte.inc'
            }
            rm(tmpInc)
        }
        ## Adjust syntax for LP solver.
        if (solver %in% c("gurobi", "lpsolveapi")) {
            env$lpobj$sense <- c("<=", env$lpobj$sense)
        }
        if (solver == "cplexapi") {
            env$lpobj$sense <- c("L", env$lpobj$sense)
        }
    } else {
        env$lpobj$rhs <- env$lpobj$rhs[-1]
        env$lpobj$sense <- env$lpobj$sense[-1]
        env$lpobj$A <- env$lpobj$A[-1, ]
    }
}

#' Minimizing violation of observational equivalence
#'
#' Given a set of IV-like estimates and the set of matrices/vectors
#' defining an LP problem, this function minimizes the violation of
#' observational equivalence under the L1 norm. The LP model must be
#' passed as an environment variable, under the entry \code{$lpobj}.
#' See \code{\link{lpSetup}}.
#' @param env environment containing the matrices defining the LP
#'     problem.
#' @param sset A list of IV-like estimates and the corresponding gamma
#'     terms.
#' @param solver string, name of the package used to solve the LP
#'     problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @param debug boolean, indicates whether or not the function should
#'     provide output when obtaining bounds. The option is only
#'     applied when \code{solver = 'gurobi'}. The output provided is
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
#' formula0 = ~ 1 + u
#' formula1 = ~ 1 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                           data = dtm,
#'                           uname = u,
#'                           as.function = FALSE)
#' polynomials1 <- polyparse(formula = formula1,
#'                           data = dtm,
#'                           uname = u,
#'                            as.function = FALSE)
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
#' ## Only the entry $sset is required
#' sSet <- sSet$sset
#'
#' ## Define additional upper- and lower-bound constraints for the LP
#' ## problem.  The code below imposes a lower bound of 0.2 and upper
#' ## bound of 0.8 on the MTRs.
#' A <- matrix(0, nrow = 22, ncol = 4)
#' A <- cbind(A, rbind(cbind(1, seq(0, 1, 0.1)),
#'                     matrix(0, nrow = 11, ncol = 2)))
#' A <- cbind(A, rbind(matrix(0, nrow = 11, ncol = 2),
#'                     cbind(1, seq(0, 1, 0.1))))
#' sense <- c(rep(">", 11), rep("<", 11))
#' rhs <- c(rep(0.2, 11), rep(0.8, 11))
#'
#' ## Construct LP object to be interpreted and solved by
#' ## lpSolveAPI. Note that an environment has to be created for the LP
#' ## object. The matrices defining the shape restrictions must be stored
#' ## as a list under the entry \code{$mbobj} in the environment.
#' lpEnv <- new.env()
#' lpEnv$mbobj <- list(mbA = A,
#'                     mbs = sense,
#'                     mbrhs = rhs)
#' ## Convert the matrices defining the shape constraints into a format
#' ## that is suitable for the LP solver.
#' lpSetup(env = lpEnv,
#'         sset = sSet,
#'         solver = "lpsolveapi")
#' ## Setup LP model so that it will minimize the criterion
#' lpSetupCriterion(env = lpEnv,
#'                 sset = sSet)
#' ## Declare any LP solver options as a list.
#' lpOptions <- optionsLpSolveAPI(list(epslevel = "tight"))
#' ## Minimize the criterion.
#' obseqMin <- criterionMin(env = lpEnv,
#'                          sset = sSet,
#'                          solver = "lpsolveapi",
#'                          solver.options = lpOptions)
#' obseqMin
#' cat("The minimum criterion is",  obseqMin$obj, "\n")
#'
#' @export
criterionMin <- function(env, sset, solver, solver.options, rescale,
                         debug = FALSE) {
    solver <- tolower(solver)
    if (solver == "gurobi") {
        if (debug && solver.options$outputflag == 1) {
            cat("\nMinimum criterion optimization statistics:\n")
            cat("------------------------------------------\n")
        }
        env$lpobj$modelsense <- "min"
        if (debug) {
            gurobi::gurobi_write(env$lpobj, "lpCriterion.mps")
            model <- env$lpobj
            save(model, file = "lpCriterion.Rdata")
            rm(model)
        }
        result   <- runGurobi(env$lpobj, solver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
        if (debug) cat("\n")
    } else if (solver == "cplexapi") {
        result   <- runCplexAPI(env$lpobj, cplexAPI::CPX_MIN, solver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
    } else if (solver == "lpsolveapi") {
        result   <- runLpSolveAPI(env$lpobj, 'min', solver.options)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
    } else if (solver == 'rmosek') {
        if (debug && solver.options$verbose == 10) {
            cat("\nMinimum criterion optimization statistics:\n")
            cat("------------------------------------------\n")
        }
        result   <- runMosek(env$lpobj, 'min', solver.options, debug)
        obseqmin <- result$objval
        optx     <- result$optx
        status   <- result$status
        if (debug) cat("\n")
    } else {
        stop(gsub('\\s+', ' ',
                  "Invalid LP solver. Option 'solver' must be either 'gurobi',
                  'cplexapi', 'rmosek', or 'lpsolveapi'."))
    }
    ## Provide nicer output
    g0sol <- optx[(2 * env$lpobj$sn + 1) :
                  (2 * env$lpobj$sn + env$lpobj$gn0)]
    g1sol <- optx[(2 * env$lpobj$sn + env$lpobj$gn0 + 1) :
                  (2 * env$lpobj$sn + env$lpobj$gn0 + env$lpobj$gn1)]
    names(g0sol) <- names(sset$gstar$g0)
    names(g1sol) <- names(sset$gstar$g1)
    if (rescale) {
        g0sol <- g0sol / env$colNorms[1:ncol(sset$s1$g0)]
        g1sol <- g1sol / env$colNorms[(ncol(sset$s1$g0) + 1):
                                      (ncol(sset$s1$g0) + ncol(sset$s1$g1))]
    }
    output <- list(obj = obseqmin,
                   x = optx,
                   g0 = g0sol,
                   g1 = g1sol,
                   status = status)
    return(output)
}

#' Obtaining TE bounds
#'
#' This function estimates the bounds on the target treatment
#' effect. The LP model must be passed as an environment variable,
#' under the entry \code{$lpobj}. See \code{\link{lpSetup}}.
#' @param env environment containing the matrices defining the LP
#'     problem.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set. This
#'     object is only used to determine the names of terms. If it is
#'     no submitted, then no names are provided to the solution
#'     vector.
#' @param noisy boolean, set to \code{TRUE} if optimization results
#'     should be displayed.
#' @param solver string, name of the package used to solve the LP
#'     problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @param smallreturnlist boolean, set to \code{TRUE} if the LP model
#'     should not be returned.
#' @param rescale boolean, set to \code{TRUE} if the MTR components
#'     should be rescaled to improve stability in the LP/QP/QCP
#'     optimization.
#' @param debug boolean, indicates whether or not the function should
#'     provide output when obtaining bounds. The option is only
#'     applied when \code{solver = 'gurobi'}. The output provided is
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
#' formula0 = ~ 1 + u
#' formula1 = ~ 1 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                           data = dtm,
#'                           uname = u,
#'                           as.function = FALSE)
#' polynomials1 <- polyparse(formula = formula1,
#'                           data = dtm,
#'                           uname = u,
#'                            as.function = FALSE)
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
#' ## Only the entry $sset is required
#' sSet <- sSet$sset
#'
#' ## Define additional upper- and lower-bound constraints for the LP
#' ## problem
#' A <- matrix(0, nrow = 22, ncol = 4)
#' A <- cbind(A, rbind(cbind(1, seq(0, 1, 0.1)),
#'                     matrix(0, nrow = 11, ncol = 2)))
#' A <- cbind(A, rbind(matrix(0, nrow = 11, ncol = 2),
#'                     cbind(1, seq(0, 1, 0.1))))
#' sense <- c(rep(">", 11), rep("<", 11))
#' rhs <- c(rep(0.2, 11), rep(0.8, 11))
#'
#' ## Construct LP object to be interpreted and solved by
#' ## lpSolveAPI. Note that an environment has to be created for the LP
#' ## object. The matrices defining the shape restrictions must be stored
#' ## as a list under the entry \code{$mbobj} in the environment.
#' lpEnv <- new.env()
#' lpEnv$mbobj <- list(mbA = A,
#'                     mbs = sense,
#'                     mbrhs = rhs)
#' ## Convert the matrices defining the shape constraints into a format
#' ## that is suitable for the LP solver.
#' lpSetup(env = lpEnv,
#'         sset = sSet,
#'         solver = "lpsolveapi")
#' ## Setup LP model so that it is solving for the bounds.
#' lpSetupBound(env = lpEnv,
#'              g0 = targetGamma$gstar0,
#'              g1 = targetGamma$gstar1,
#'              sset = sSet,
#'              criterion.tol = 0,
#'              criterion.min = 0,
#'              solver = "lpsolveapi")
#' ## Declare any LP solver options as a list.
#' lpOptions <- optionsLpSolveAPI(list(epslevel = "tight"))
#' ## Obtain the bounds.
#' bounds <- bound(env = lpEnv,
#'                 sset = sSet,
#'                 solver = "lpsolveapi",
#'                 solver.options = lpOptions)
#' cat("The bounds are [",  bounds$min, ",", bounds$max, "].\n")
#'
#' @export
bound <- function(env, sset, solver,
                  solver.options, noisy = FALSE,
                  smallreturnlist = FALSE,
                  rescale = FALSE,
                  debug = FALSE) {
    solver <- tolower(solver)
    ## Obtain lower and upper bounds
    if (solver == "gurobi") {
        if (debug && solver.options$outputflag == 1) {
            cat("\nLower bound optimization statistics:\n")
            cat("------------------------------------\n")
        }
        if (debug == TRUE){
            gurobi::gurobi_write(env$lpobj, "lpBound.mps")
            model <- env$lpobj
            save(model, file = "lpBound.Rdata")
            rm(model)
        }
        env$lpobj$modelsense <- "min"
        minresult <- runGurobi(env$lpobj, solver.options)
        min <- minresult$objval
        minstatus <- minresult$status
        minoptx <- minresult$optx
        if (debug && solver.options$outputflag == 1) {
            cat("\nUpper bound optimization statistics:\n")
            cat("------------------------------------\n")
        }
        env$lpobj$modelsense <- "max"
        maxresult <- runGurobi(env$lpobj, solver.options)
        max <- maxresult$objval
        maxstatus <- maxresult$status
        maxoptx <- maxresult$optx
        if (debug) cat("\n")
    } else if (solver == "cplexapi") {
        minresult <- runCplexAPI(env$lpobj, cplexAPI::CPX_MIN, solver.options)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runCplexAPI(env$lpobj, cplexAPI::CPX_MAX, solver.options)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
    } else if (solver == "lpsolveapi") {
        minresult <- runLpSolveAPI(env$lpobj, 'min', solver.options)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        maxresult <- runLpSolveAPI(env$lpobj, 'max', solver.options)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
    } else if (solver == 'rmosek') {
        if (debug && solver.options$verbose == 10) {
            cat("\nLower bound optimization statistics:\n")
            cat("------------------------------------\n")
        }
        minresult <- runMosek(env$lpobj, 'min', solver.options, debug)
        min       <- minresult$objval
        minoptx   <- minresult$optx
        minstatus <- minresult$status
        if (debug && solver.options$verbose == 10) {
            cat("\nUpper bound optimization statistics:\n")
            cat("------------------------------------\n")
        }
        maxresult <- runMosek(env$lpobj, 'max', solver.options)
        max       <- maxresult$objval
        maxoptx   <- maxresult$optx
        maxstatus <- maxresult$status
        if (debug) cat("\n")
    } else {
        stop(gsub('\\s+', ' ',
                  "Invalid LP solver. Option 'solver' must be either 'gurobi',
                  'cplexapi', 'rmosek', or 'lpsolveapi'."))
    }
    env$lpobj$modelsense <- NULL
    ## Return error codes, if any
    if (maxstatus %in% c(2, 3, 4, 5, 9) || minstatus %in% c(2, 3, 4, 5, 9)) {
        return(list(error = TRUE,
                    max = max,
                    maxstatus = maxstatus,
                    min = min,
                    minstatus = minstatus))
    }
    ming0 <- minoptx[(2 * env$lpobj$sn + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0)]
    ming1 <- minoptx[(2 * env$lpobj$sn + env$lpobj$gn0 + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0 + env$lpobj$gn1)]
    maxg0 <- maxoptx[(2 * env$lpobj$sn + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0)]
    maxg1 <- maxoptx[(2 * env$lpobj$sn + env$lpobj$gn0 + 1) :
                     (2 * env$lpobj$sn + env$lpobj$gn0 + env$lpobj$gn1)]
    if (rescale) {
        ming0 <- ming0 / env$colNorms[1:ncol(sset$s1$g0)]
        ming1 <- ming1 / env$colNorms[(ncol(sset$s1$g0) + 1):
                                      (ncol(sset$s1$g0) + ncol(sset$s1$g1))]
        maxg0 <- maxg0 / env$colNorms[1:ncol(sset$s1$g0)]
        maxg1 <- maxg1 / env$colNorms[(ncol(sset$s1$g0) + 1):
                                      (ncol(sset$s1$g0) + ncol(sset$s1$g1))]
    }
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
    ## Name the coefficients
    if (is.matrix(sset$s1$g0)) {
        names(ming0) <- names(maxg0) <- colnames(sset$s1$g0)
        names(ming1) <- names(maxg1) <- colnames(sset$s1$g1)
    } else {
        names(ming0) <- names(maxg0) <- names(sset$s1$g0)
        names(ming1) <- names(maxg1) <- names(sset$s1$g1)
    }
    ## Return output
    output <- list(max = max,
                   maxg0 = maxg0,
                   maxg1 = maxg1,
                   maxresult = maxresult,
                   maxstatus = maxstatus,
                   min = min,
                   ming0 = ming0,
                   ming1 = ming1,
                   minresult = minresult,
                   minstatus = minstatus,
                   error = FALSE)
    if (rescale) {
        output$norms <- env$colNorms
    }
    if (!smallreturnlist) {
        if (solver != 'rmosek') output$model = env$lpobj
        if (solver == 'rmosek') {
            output$model <- minresult$prob
            output$model$sense <- NULL
        }
    }
    return(output)
}


#' Running Gurobi LP solver
#'
#' This function solves the LP problem using the Gurobi package. The
#' object generated by \code{\link{lpSetup}} is compatible with the
#' \code{gurobi} function. See \code{\link{runCplexAPI}} for
#' additional error code labels.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from Gurobi. This includes the
#'     objective value, the solution vector, and the optimization
#'     status (status of \code{1} indicates successful optimization) .
runGurobi <- function(lpobj, solver.options) {
    result <- gurobi::gurobi(lpobj, solver.options)
    status <- 0
    if (result$status == "OPTIMAL") status <- 1
    if (result$status == "INFEASIBLE") status <- 2
    if (result$status == "INF_OR_UNBD") status <- 3
    if (result$status == "UNBOUNDED") status <- 4
    if (result$status == "NUMERIC") status <- 5
    if (result$status == "SUBOPTIMAL") status <- 6
    optx <- result$x
    return(list(objval = result$objval,
                optx = result$x,
                status = status))
}


#' Running cplexAPI LP solver
#'
#' This function solves the LP problem using the cplexAPI package. The
#' object generated by \code{\link{lpSetup}} is not compatible with
#' the \code{cplexAPI} functions. This function adapts the object to
#' solve the LP problem. See \code{\link{runGurobi}} for additional
#' error code labels.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param lpdir input either CPX_MAX or CPX_MIN, which sets the LP
#'     problem as a maximization or minimization problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from CPLEX. This includes the
#'     objective value, the solution vector, and the optimization
#'     status (status of \code{1} indicates successful optimization).
runCplexAPI <- function(lpobj, lpdir, solver.options) {
    ## Declare environment and set options
    env  <- cplexAPI::openEnvCPLEX()
    prob <- cplexAPI::initProbCPLEX(env)
    cplexAPI::chgProbNameCPLEX(env, prob, "sample")
    if (!is.null(solver.options)) {
        for(i in seq(length(solver.options))) {
            eval(parse(text = solver.options[[i]]))
        }
    }
    ## Declare LP prblem
    sense <- lpobj$sense
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
                                lb = lpobj$lb,
                                ub = lpobj$ub)
    cplexAPI::lpoptCPLEX(env, prob)
    solution <- cplexAPI::solutionCPLEX(env, prob)
    cplexAPI::delProbCPLEX(env, prob)
    cplexAPI::closeEnvCPLEX(env)
    status <- 0
    if (typeof(solution) == "S4") {
        if (attr(solution, "class") == "cplexError") {
            status <- 5
            solution <- list()
            solution$objval <- NA
            solution$x <- NA
        }
    }  else {
        if (solution$lpstat == 1) status <- 1
        if (solution$lpstat == 2) status <- 4
        if (solution$lpstat == 3) status <- 2
        if (solution$lpstat == 4) status <- 3
        if (solution$lpstat == 5) status <- 7
        if (solution$lpstat == 6) status <- 6
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
#' adapts the object to solve the LP problem. See
#' \code{\link{runGurobi}} and \code{\link{runCplexAPI}} for
#' additional error code labels.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param modelsense input either 'max' or 'min' which sets the LP
#'     problem as a maximization or minimization problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from \code{lpSolveAPI}. This includes
#'     the objective value, the solution vector, and the optimization
#'     status (status of \code{1} indicates successful optimization).
runLpSolveAPI <- function(lpobj, modelsense, solver.options) {
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
    if (!is.null(solver.options)) {
        eval(solver.options)
    }
    lpSolveAPI::set.bounds(lprec = lpmodel,
                           lower = lpobj$lb,
                           upper = lpobj$ub)
    solved <- lpSolveAPI::solve.lpExtPtr(lpmodel)
    status <- 0
    if (solved == 0) status <- 1
    if (solved == 1) status <- 6
    if (solved == 2) status <- 2
    if (solved == 3) status <- 4
    if (solved == 5) status <- 5
    return(list(objval = lpSolveAPI::get.objective(lpmodel),
                optx   = lpSolveAPI::get.variables(lpmodel),
                status = status))
}

#' Running Rmosek
#'
#' This function solves the LP problem using the \code{Rmosek}
#' package. The object generated by \code{\link{lpSetup}} is not
#' compatible with the \code{Rmosek} functions. This function
#' adapts the object to solve the LP problem. See
#' \code{\link{runGurobi}} and \code{\link{runCplexAPI}} for
#' additional error code labels.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param modelsense input either 'max' or 'min' which sets the LP
#'     problem as a maximization or minimization problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from \code{lpSolveAPI}. This includes
#'     the objective value, the solution vector, and the optimization
#'     status (status of \code{1} indicates successful optimization).
runMosek <- function(lpobj, modelsense, solver.options, debug = FALSE) {
    qcp <- !is.null(lpobj$Q)
    qcqp <- !is.null(lpobj$quadcon)
    if (qcqp) {
        ## Construct the transformation matrix
        Qc <- lpobj$quadcon[[1]]$Qc * 2
        QdList <- eigen(Qc, symmetric = TRUE)
        QdList$values <- sapply(QdList$values, function(x) max(x, 0))
        Qd <- diag(sqrt(QdList$values)) %*% t(QdList$vectors)
        ## Update the old linear constraints to include the new variables
        nvars <- length(lpobj$obj)
        QSlack1 <- c(rep(0, times = nvars), c(1, 0), rep(0, times = nvars))
        QSlack2 <- c(lpobj$quadcon[[1]]$q, 0, 1, rep(0, times = nvars))
        lpobj$A <- rbind(cbind(lpobj$A, Matrix::Matrix(0, nrow = nrow(lpobj$A),
                                                       ncol = nvars + 2)),
                         QSlack1,
                         QSlack2,
                         cbind(Qd,
                               Matrix::Matrix(0, nrow = nvars, ncol = 2),
                               -diag(nvars)))
        lpobj$rhs <- c(lpobj$rhs, 1, lpobj$quadcon[[1]]$rhs, rep(0, nvars))
        lpobj$sense <- c(lpobj$sense, '=', '=', rep('=', nvars))
        lpobj$ub <- c(lpobj$ub, 1, Inf, rep(Inf, nvars))
        lpobj$lb <- c(lpobj$lb, 1, 0, rep(-Inf, nvars))
        ## Update objective to include new variables
        lpobj$obj <- c(lpobj$obj, rep(0, nvars + 2))
    }
    prob <- list()
    tmpSense <- lpobj$sense
    tmpBlc <- lpobj$rhs
    tmpBuc <- lpobj$rhs
    tmpBuc[which(tmpSense == '>=')] <- Inf
    tmpBlc[which(tmpSense == '<=')] <- -Inf
    prob$bc <- rbind(blc = tmpBlc,
                     buc = tmpBuc)
    prob$bx <- rbind(blx = lpobj$lb,
                     bux = lpobj$ub)
    prob$A <- lpobj$A
    prob$c <- lpobj$obj
    prob$sense <- modelsense
    if (qcp) {
        tripletQ <- matrixTriplets(lpobj$Q, lower = TRUE)
        prob$qobj$i <- tripletQ$rows
        prob$qobj$j <- tripletQ$columns
        prob$qobj$v <- tripletQ$values * 2
    }
    if (qcqp) {
        ## Declare conic constraints
        prob$cones <- matrix(list(), nrow = 2, ncol = 1)
        rownames(prob$cones) <- c('type', 'sub')
        prob$cones[, 1] <- list('RQUAD', seq((nvars + 1), (2 * nvars + 2)))
    }
    ## Include options
    if ('dparam' %in% names(solver.options)) {
        prob$dparam <- solver.options$dparam
        solver.options$dparam <- NULL
    }
    ## Export model if debugging
    if (debug == TRUE){
        if (qcp) {
            Rmosek::mosek_write(prob, "lpCriterion.mps", solver.options)
            save(prob, file = "lpCriterion.Rdata")
        }
        if (qcqp) {
            Rmosek::mosek_write(prob, "lpBound.mps", solver.options)
            save(prob, file = "lpBound.Rdata")
        }
    }
    ## Estimate
    result <- Rmosek::mosek(prob, solver.options)
    response <- result$response$code
    if (is.null(result$sol)) {
        objval <- NA
        optx <- NA
        responsecode <- result$response$msg
        status <- 0
        solutionstatus <- NA
        problemstatus <- NA
    } else {
        if (!is.null(result$sol$itr$xc)) {
            optx <- result$sol$itr$xx
            solutionstatus <- result$sol$itr$solsta
            problemstatus <- result$sol$itr$prosta
        } else {
            optx <- result$sol$bas$xx
            solutionstatus <- result$sol$bas$solsta
            problemstatus <- result$sol$bas$prosta
        }
        if (solutionstatus == 'OPTIMAL') {
            status <- 1
        } else if (solutionstatus == 'DUAL_INFEASIBLE_CER') {
            status <- 4
        } else if (solutionstatus == 'PRIMAL_INFEASIBLE_CER') {
            status <- 2
        } else if (solutionstatus == 'UNKNOWN') {
            if (response == 0) {
                status <- 8
            } else if (response == 10000) {
                status <- 10
            } else if (response == 10006) {
                status <- 6
            } else if (response == 10025) {
                status <- 5
            } else {
                status <- 0
            }
        }
        if (!is.null(optx)) {
            if (!qcp) {
                objval <- sum(optx * lpobj$obj)
            } else {
                objval <- t(optx) %*% lpobj$Q %*% optx +
                    t(optx) %*% lpobj$obj
            }
            if (qcqp) {
                optx <- optx[1:nvars]
            }
        }
    }
    return(list(objval = as.numeric(objval),
                optx = as.vector(optx),
                response = response,
                status = status,
                solutionstatus = solutionstatus,
                problemstatus = problemstatus,
                prob = prob))
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
#' \code{solver} is set to \code{Gurobi}. This function really
#' implements some default values, and accounts for the \code{debug}
#' option.
#' @param options list. The list should be structured the same way as
#'     if one were using the \code{gurobi} library directly. That is,
#'     the name of each item must be the name of the option, and is
#'     case sensitive. The value assigned to each item is the value to
#'     set the option to.
#' @param debug boolean, indicates whether or not the function should
#'     provide output when obtaining bounds. The output provided is
#'     the same as what the Gurobi API would send to the console.
#' @return list, the set of options declared by the user, including
#'     some additional default values (if not assigned by the user)
#'     and accounting for \code{debug}.
#' @export
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
        options$presolve <- 1
    }
    return(options)
}

#' Function to parse options for Gurobi
#'
#' This function constructs a list of options to be parsed when
#' \code{solver} is set to \code{Rmosek}. This function really
#' implements the default feasibility tolerances.
#' @param options list. Each set of options should be passed as a
#'     list, with the name of each entry being the name of the class
#'     of options. For example, options for double parameters should
#'     be contained in the entry\code{dparam = list(BASIS_TOL_X = 1e-06)}.
#' @param debug boolean, indicates whether or not the function should
#'     provide output when obtaining bounds. The output provided is
#'     the same as what Mosek would send to the console.
#' @return list, the set of options declared by the user, including
#'     some additional default values.
#' @export
optionsRmosek <- function(options, debug) {
    if (is.null(options$dparam$ANA_SOL_INFEAS_TOL)) {
        options$dparam$ANA_SOL_INFEAS_TOL <- 1e-06
    }
    if (is.null(options$dparam$BASIS_TOL_X)) {
        options$dparam$BASIS_TOL_X <- 1e-06
    }
    if (is.null(options$verbose)) {
        if (debug) options$verbose <- 10
        if (!debug) options$verbose <- 0
    }
    return(options)
}

#' Function to parse options for lp_solve
#'
#' This function constructs a list of options to be parsed when
#' \code{solver} is set to \code{lpsolveapi}. The options permitted
#' are those that can be set via \code{lpSolveAPI::lp.control}, and
#' should be passed as a named list (e.g. \code{list(epslevel =
#' "tight")}).
#' @param options list. The name of each item must be the name of the
#'     option, and is case sensitive. The value assigned to each item
#'     is the value to set the option to. The \code{lprec} argument
#'     should always be omitted.
#' @return string, the command to be evaluated to implement the
#'     options.
#' @export
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
#' \code{solver} is set to \code{cplexapi}.
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
#'     passed as the parameter value (e.g. \code{list(setDefaultParm
#'     = NA)}).
#' @return list, each element being the command to evaluate to
#'     implement an option.
#' @export
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
#' This function constructs a string to be parsed when \code{solver}
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

#' Constructing QCQP problem
#'
#' This function is only used when the direct MTR regression procedure
#' is used. This function simply constructs the quadratic constraint,
#' and adds it to the LP problem defined by the linear optimization problem
#' for the bounds and the linear shape constraints.
#'
#' @param env environment containing the matrices defining the LP
#'     problem.
#' @param sset A list containing the covariats and outcome variable
#'     for the direct MTR regression.
#' @param rescale boolean, set to \code{TRUE} if the MTR components
#'     should be rescaled to improve stability in the LP/QP/QCP
#'     optimization.
#' @export
qpSetup <- function(env, sset, rescale = TRUE) {
    ## Construct the constraint vectors and matrices
    drY <- sset$s1$ys
    drX <- cbind(sset$s1$g0, sset$s1$g1)
    drN <- length(drY)
    if (rescale) {
        drX <- sweep(x = drX, MARGIN = 2, STATS = env$colNorms, FUN = '/')
        normY <- sqrt(sum(drY^2))
        drN <- drN * normY
    }
    qc <- list()
    qc$q <- as.vector(-2 * t(drX) %*% drY) / drN
    qc$Qc <- Matrix::Matrix(t(drX) %*% drX / drN, sparse = TRUE)
    qc$sense <- '<'
    qc$name <- 'SSR'
    ## Store the quadratic component
    env$quad <- qc
    ## Store the SSY
    env$ssy <- sum(drY^2)
    ## Store the regression matrices
    env$drY <- drY
    env$drX <- drX
    env$drN <- drN
    if (rescale) env$normY <- normY
}

#' Configure QCQP problem to find minimum criterion
#'
#' This function sets up the objective function for minimizing the
#' criterion. The QCQP model must be passed as an environment variable,
#' under the entry \code{$lpobj}. See \code{\link{qpSetup}}.
#' @param env The LP environment
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
qpSetupCriterion <- function(env) {
    env$lpobj$obj <- env$quad$q
    env$lpobj$Q <- env$quad$Qc
    env$lpobj$modelsense <- 'min'
    env$lpobj$ub <- rep(Inf, times = length(env$quad$q))
    env$lpobj$lb <- rep(-Inf, times = length(env$quad$q))
}

#' Constructing QCQP problem for bounding
#'
#' This function is only used when the direct MTR regression procedure
#' is used. This function simply constructs the quadratic constraint,
#' and adds it to the LP problem defined by the linear optimization problem
#' for the bounds and the linear shape constraints.
#'
#' @param env environment containing the matrices defining the LP
#'     problem.
#' @param g0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param g1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param criterion.tol non-negative scalar, determines how much the
#'     quadratic constraint should be relaxed by. If set to 0, the
#'     constraint is not relaxed at all.
#' @param criterion.min minimum of (SSR - SSY) of a linear regression
#'     with shape constraints.
#' @param rescale boolean, set to \code{TRUE} if the MTR components
#'     should be rescaled to improve stability in the LP/QP/QCP
#'     optimization.
#' @param setup boolean, set to \code{TRUE} if the QP problem should
#'     be set up for solving the bounds, which includes the quadratic
#'     constraint. Set to \code{FALSE} if the quadratic constraint
#'     should be removed.
#' @return A list of matrices and vectors necessary to define an LP
#'     problem for Gurobi.
#' @export
qpSetupBound <- function(env, g0, g1,
                         criterion.tol,
                         criterion.min,
                         rescale = FALSE,
                         setup = TRUE) {
    if (setup) {
        ## Prepare objective
        env$lpobj$obj <- c(g0, g1)
        if (rescale) {
            env$lpobj$obj <- env$lpobj$obj / env$colNorms
        }
        env$lpobj$Q <- NULL
        ## Add in the quadratic constraint, accounting for how
        ## criterion.min excludes the SSY.
        env$quad$rhs <- (env$ssy / env$drN) * criterion.tol +
            criterion.min * (1 + criterion.tol)
        env$lpobj$quadcon <- list(env$quad)
    } else {
        env$lpobj$obj <- NULL
        env$quad$rhs <- NULL
        env$lpobj$quadcon <- NULL
        env$lpobj$ub <- NULL
        env$lpobj$lb <- NULL
    }
}

#' Configure QP environment for diagnostics
#'
#' This function separates the shape constraints from the QP
#' environment. That way, the model can be solved without any shape
#' constraints, which is the primary cause of infeasibility. This is
#' done in order to check which shape constraints are causing the
#' model to be infeasible. The QP model must be passed as an
#' environment variable, under the entry \code{$lpobj}. See
#' \code{\link{lpSetup}}.
#' @param env The LP environment
#' @param rescale boolean, set to \code{TRUE} if the MTR components
#'     should be rescaled to improve stability in the LP/QP/QCP
#'     optimization.
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
qpSetupInfeasible <- function(env, rescale) {
    ## Separate shape constraint objects
    env$mbobj$mbA <- env$lpobj$A
    env$mbobj$mbrhs <- env$lpobj$rhs
    env$mbobj$mbsense <- env$lpobj$rhs
    ## Reduce lpobjects so they do not contain any shape constraints
    env$lpobj$A <- matrix(0, nrow = 1, ncol = ncol(env$mbobj$mbA))
    env$lpobj$rhs <- 0
    env$lpobj$sense <- '='
}

#' Convert matrix into triplet form
#'
#' This function converts matrices into triplet form for Mosek.  This
#' is required in order to declare quadratic programming problems and
#' second-order cone programming problems.
#'
#' @param mat A matrix.
#' @param lower Boolean, set to \code{TRUE} if matrix is symmetric,
#'     and only its lower triangle should be returned.
#' @return A list containing vectors of row and column indexes, and
#'     matrix values.
matrixTriplets <- function(mat, lower = TRUE) {
    dim <- dim(mat)
    sparseList <- apply(mat, 2, function(x) {
        pos <- which(x != 0)
        vals <- x[x != 0]
        return(list(pos = pos,
                    vals = vals))
    })
    ## If lower = TRUE, that means return entries for the lower
    ## triangle only.
    if (lower) {
        sparseList <- lapply(seq(dim[2]), function(x) {
            pos <- sparseList[[x]]$pos
            vals <- sparseList[[x]]$vals
            if (length(pos) > 0) {
                vals <- vals[pos >= x]
                pos <- pos[pos >= x]

            }
            return(list(pos = pos,
                        vals = vals))
        })
    }
    columnIndex <- Reduce(c, lapply(seq(dim[2]),
                                    FUN = function(x) {
                                        n <- length(sparseList[[x]]$pos)
                                        if (n > 0) return(rep(x, times = n))
                                    }))
    rowIndex <- Reduce(c, lapply(sparseList, function(x) x$pos))
    values <- Reduce(c, lapply(sparseList, function(x) x$vals))
    return(list(columns = columnIndex,
                rows = rowIndex,
                values = values,
                dim = dim,
                lower = lower))
}
