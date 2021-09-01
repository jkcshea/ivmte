#' Audit procedure
#'
#' This is the wrapper for running the entire audit procedure. This
#' function sets up the LP/QCQP problem of minimizing criterion.  for
#' the set of IV-like estimands, while satisfying boundedness and
#' monotonicity constraints declared by the user. Rather than enforce
#' that boundedness and monotonicity hold across the entire support of
#' covariates and unobservables, this procedure enforces the
#' conditions over a grid of points. This grid corresponds to the set
#' of values the covariates can take, and a set of values of the
#' unobservable term. The size of this grid is specified by the user
#' in the function arguments. The procedure first estimates the bounds
#' while imposing the shape constraints for an initial subset of
#' points in the grid. The procedure then goes on to check ('audit')
#' whether the constraints are satisfied over the entire grid. Any
#' point where either the boundedness or monotonicity constraints are
#' violated are incorporated into the initial grid, and the process is
#' repeated until the audit no longer finds any violations, or until
#' some maximum number of iterations is reached.
#'
#' @param audit.grid list, contains the \code{A} matrix used in the
#'     audit for the original sample, as well as the RHS vector used
#'     in the audit from the original sample.
#' @param pm0 A list of the monomials in the MTR for the control
#'     group.
#' @param pm1 A list of the monomials in the MTR for the treated
#'     group.
#' @param m1.ub.default boolean, default set to FALSE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m0.ub.default boolean, default set to FALSE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param mte.ub.default boolean, default set to FALSE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m1.lb.default boolean, default set to FALSE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m0.lb.default boolean, default set to FALSE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param mte.lb.default boolean, default set to FALSE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param equal.coef0 character, a vector containing all the terms in
#'     \code{m0} that should have the same coefficients in
#'     \code{m1}. The order of the variables must match those of
#'     \code{equal.coef1}, which contains all the corresponding terms
#'     in \code{m1}. The reason the terms are entered separately for
#'     \code{m0} and \code{m1} is because the spline terms may be
#'     named differently across treatment and control groups.
#' @param equal.coef1 character, a vector containing all the terms in
#'     \code{m1} that should have the same coefficients in
#'     \code{m0}. See the description for \code{equal.coef0} for more
#'     details.
#' @param sset a list containing the point estimates and gamma moments
#'     for each IV-like specification.
#' @param gstar0 set of expectations for each terms of the MTR for the
#'     control group, corresponding to the target parameter.
#' @param gstar1 set of expectations for each terms of the MTR for the
#'     control group, corresponding to the target parameter.
#' @param orig.sset list, only used for bootstraps. The list contains
#'     the gamma moments for each element in the S-set, as well as the
#'     IV-like coefficients.
#' @param orig.criterion numeric, only used for bootstraps. The scalar
#'     corresponds to the minimum observational equivalence criterion
#'     from the original sample.
#' @param rescale boolean, set to \code{TRUE} if the MTR components
#'     should be rescaled to improve stability in the LP/QCQP
#'     optimization.
#'
#' @inheritParams ivmteEstimate
#'
#' @return a list. Included in the list are estimates of the treatment
#'     effect bounds; the minimum violation of observational
#'     equivalence of the set of IV-like estimands; the list of
#'     matrices and vectors defining the LP/QCQP problem; the points used
#'     to generate the audit grid, and the points where the shape
#'     constraints were violated.
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
#' ## obtained from this object by the function 'genSSet'
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## If splines are interacted with other variables, the
#' ## 'interactSplines' should be used.
#' ## splinesList <- interactSplines(splinesobj = splinesList,
#' ##                               m0 = formula0,
#' ##                               m1 = formula1,
#' ##                               data = data,
#' ##                               uname = 'u')
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' polynomials1 <- polyparse(formula = formula1,
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
#' ## Perform audit procedure and return bounds
#' audit(data = dtm,
#'       uname = u,
#'       m0 = formula0,
#'       m1 = formula1,
#'       pm0 = polynomials0,
#'       pm1 = polynomials1,
#'       splinesobj = splinesList,
#'       vars_data = colnames(dtm),
#'       vars_mtr = "u",
#'       terms_mtr0 = "u",
#'       terms_mtr1 = "u",
#'       sset = sSet$sset,
#'       gstar0 = targetGamma$gstar0,
#'       gstar1 = targetGamma$gstar1,
#'       m0.inc = TRUE,
#'       m1.dec = TRUE,
#'       m0.lb = 0.2,
#'       m1.ub = 0.8,
#'       audit.max = 5,
#'       solver = "lpSolveAPI")
#'
#' @export
audit <- function(data, uname, m0, m1, pm0, pm1, splinesobj,
                  vars_mtr, terms_mtr0, terms_mtr1, vars_data,
                  initgrid.nu = 20, initgrid.nx = 20,
                  audit.nx = 2500, audit.nu = 25, audit.add = 100,
                  audit.max = 25, audit.tol,
                  audit.grid = NULL,
                  m1.ub, m0.ub, m1.lb, m0.lb, mte.ub, mte.lb,
                  m1.ub.default = FALSE,
                  m0.ub.default = FALSE,
                  mte.ub.default = FALSE,
                  m1.lb.default = FALSE,
                  m0.lb.default = FALSE,
                  mte.lb.default = FALSE,
                  m0.dec = FALSE, m0.inc = FALSE,
                  m1.dec = FALSE, m1.inc = FALSE,
                  mte.dec = FALSE, mte.inc = FALSE,
                  equal.coef0, equal.coef1,
                  sset, gstar0, gstar1,
                  orig.sset = NULL, orig.criterion = NULL,
                  criterion.tol = 1e-4,
                  solver, solver.options, solver.presolve,
                  solver.options.criterion, solver.options.bounds,
                  rescale = TRUE,
                  smallreturnlist = FALSE,
                  noisy = TRUE, debug = FALSE) {
    call  <- match.call()
    solver <- tolower(solver)
    ## Determine if whether IV-like moments or direct MTR regression
    ## will be used.
    direct <- (length(sset) == 1 & 'SSR' %in% names(sset$s1))
    if (!direct) rescale <- FALSE
    ## Set the audit tolerance
    if (!hasArg(audit.tol)) {
        if (solver == "gurobi") {
            if (hasArg(solver.options) &&
                !is.null(solver.options$FeasibilityTol)) {
                audit.tol <- solver.options$FeasibilityTol
            } else if (hasArg(solver.options.bounds) &&
                       !is.null(solver.options.bounds$FeasibilityTol)) {
                audit.tol <- solver.options.bounds$FeasibilityTol
            } else {
                audit.tol <- 1e-06
            }
        }
        if (solver == "cplexapi") {
            if (hasArg(solver.options)) {
                audit.tol <- optionsCplexAPITol(solver.options)
            } else if (hasArg(solver.options.bounds)){
                audit.tol <- optionsCplexAPITol(solver.options.bounds)
            } else {
                audit.tol <- 1e-06
            }
        }
        if (solver == 'rmosek') {
            if (hasArg(solver.options) &&
                !is.null(solver.options$dparam$ANA_SOL_INFEAS_TOL)) {
                audit.tol <- solver.options$dparam$ANA_SOL_INFEAS_TOL
            } else if (hasArg(solver.options.bounds) &&
                       !is.null(solver.options.bounds$solver.options$dparam$ANA_SOL_INFEAS_TOL)) {
                audit.tol <- solver.options.bounds$solver.options$dparam$ANA_SOL_INFEAS_TOL
            } else {
                audit.tol <- 1e-06
            }
        }
        if (solver == "lpsolveapi") {
            audit.tol <- 1e-06
        }
    }
    ## Organize solver options
    if (solver == "gurobi") {
        ## Construct default options
        solver.options.default <- list(dualreductions = 1,
                                         FeasibilityTol = 1e-06)
        if (debug)  solver.options.default$outputflag <- 1
        if (!debug) solver.options.default$outputflag <- 0
        if (hasArg(solver.presolve)) solver.options.default$presolve <-
                                           as.integer(solver.presolve)
        ## Prepare user options
        if (hasArg(solver.options)) {
            solver.options <- optionsGurobi(solver.options, debug)
            if (hasArg(solver.presolve)) solver.options$presolve <-
                                               as.integer(solver.presolve)
            solver.options.criterion <- solver.options
            solver.options.bounds <- solver.options
        } else {
            if (hasArg(solver.options.criterion)) {
                solver.options.criterion <-
                    optionsGurobi(solver.options.criterion, debug)
                if (hasArg(solver.presolve)) {
                    solver.options.criterion$presolve <-
                        as.integer(solver.presolve)
                }
            } else {
                solver.options.criterion <- solver.options.default
            }
            if (hasArg(solver.options.bounds)) {
                solver.options.bounds <-
                    optionsGurobi(solver.options.bounds, debug)
                if (hasArg(solver.presolve)) {
                    solver.options.bounds$presolve <-
                        as.integer(solver.presolve)
                }
            } else {
                solver.options.bounds <- solver.options.default
            }
        }
        ## Turn off non-convex option if using direct regression
        if (direct) {
            solver.options.criterion$nonconvex <- 0
            solver.options.bounds$nonconvex <- 0
        }
    } else {
        if (solver == "cplexapi") {
            if (hasArg(solver.options)) {
                solver.options.criterion <- optionsCplexAPI(solver.options)
                solver.options.bounds <- optionsCplexAPI(solver.options)
            } else {
                if (hasArg(solver.options.criterion)) {
                    solver.options.criterion <-
                        optionsCplexAPI(solver.options.criterion)
                } else {
                    solver.options.criterion <-
                        optionsCplexAPI(list(setDblParmCPLEX =
                                                 c(parm = 1016,
                                                   value = 1e-06)))
                }
                if (hasArg(solver.options.bounds)) {
                    solver.options.bounds <-
                        optionsCplexAPI(solver.options.bounds)
                } else {
                    solver.options.bounds <-
                        optionsCplexAPI(list(setDblParmCPLEX =
                                                 c(parm = 1016,
                                                   value = 1e-06)))
                }
            }
        }
        if (solver == "rmosek") {
            ## Construct default options
            solver.options.default <- list()
            solver.options.default$dparam <- list(ANA_SOL_INFEAS_TOL = 1e-06,
                                                  BASIS_TOL_X = 1e-06)
            if (debug) solver.options.default$verbose <- 10
            if (!debug) solver.options.default$verbose <- 0
            ## Prepare user options
            if (hasArg(solver.options)) {
                solver.options <- optionsRmosek(solver.options, debug)
                solver.options.criterion <- solver.options
                solver.options.bounds <- solver.options
            } else {
                if (hasArg(solver.options.criterion)) {
                    solver.options.criterion <-
                        optionsRmosek(solver.options.criterion, debug)
                } else {
                    solver.options.criterion <- solver.options.default
                }
                if (hasArg(solver.options.bounds)) {
                    solver.options.bounds <-
                        optionsRmosek(solver.options.bounds, debug)
                } else {
                    solver.options.bounds <- solver.options.default
                }
            }
        }
        if (solver == "lpsolveapi") {
            if (hasArg(solver.options)) {
                solver.options.criterion <-
                    optionsLpSolveAPI(solver.options)
                solver.options.bounds <- optionsLpSolveAPI(solver.options)
            } else {
                if (hasArg(solver.options.criterion)) {
                    solver.options.criterion <-
                        optionsLpSolveAPI(solver.options.criterion)
                } else {
                    solver.options.criterion <-
                        optionsLpSolveAPI(list(epslevel = "tight"))
                }
                if (hasArg(solver.options.bounds)) {
                    solver.options.bounds <-
                        optionsLpSolveAPI(solver.options.bounds)
                } else {
                    solver.options.bounds <-
                        optionsLpSolveAPI(list(epslevel = "tight"))
                }
            }
        }
    }
    ## Clean boolean terms
    terms_mtr0 <- parenthBoolean(terms_mtr0)
    terms_mtr1 <- parenthBoolean(terms_mtr1)
    ## Obtain name of unobservable variable
    if(hasArg(uname)) {
        ## Convert uname into a string
        uname <- deparse(substitute(uname))
        uname <- gsub("~", "", uname)
        uname <- gsub("\\\"", "", uname)
    } else {
        uname <- "u"
    }
    ## Organize variables
    monov <- uname ## `monov' is a placeholder name for the monotone
                   ## variable. I use this in case I want to
                   ## generalize the monotonciity restrictions to
                   ## other covariates
    xvars <- unique(vars_mtr)
    xvars <- xvars[xvars != uname]
    xvars <- xvars[xvars %in% vars_data]
    otherx  <- xvars[xvars != monov]
    ## Begin performing the audit
    prevbound <- c(-Inf, Inf)
    existsolution <- FALSE
    audit_count <- 1
    ## Generate the components for the audit grid, and generate the
    ## initial grid
    monoboundAlist <- c('sset', 'gstar0', 'gstar1',
                        'm1.ub', 'm0.ub',
                        'm1.lb', 'm0.lb',
                        'mte.ub', 'mte.lb',
                        'm0.dec', 'm0.inc',
                        'm1.dec', 'm1.inc',
                        'mte.dec', 'mte.inc',
                        'pm0', 'pm1')
    if (audit_count == 1) {
        if (!direct) sn <- length(sset)
        if (direct) sn <- 0
        if (is.null(audit.grid)) {
            ## Generate the underlying X grid for the audit
            if (length(xvars) == 0) {
                noX <- TRUE
                support <- NULL
            } else {
                noX <- FALSE
                support <- unique(data[, xvars])
                ## Check if support is vector or matrix; it can be a vector if
                ## there is only one X term
                if (is.null(dim(support))) {
                    support <- data.frame(support)
                    colnames(support) <- xvars
                }
                rownames(support) <- seq(1, nrow(support))
                ## Select first iteration of the grid
                full_index <- seq(1, nrow(support))
                initgrid.nx <- min(initgrid.nx, nrow(support))
                audit.nx <- min(audit.nx, nrow(support))
                a_grid_index <- sort(sample(full_index,
                                            audit.nx,
                                            replace = FALSE,
                                            prob = replicate(length(full_index),
                                            (1 / length(full_index)))))
                ## Restrict support the audit points
                support <- support[a_grid_index, ]
                if (is.null(dim(support))) {
                    support <- data.frame(support)
                    colnames(support) <- xvars
                }
                rownames(support) <- seq(nrow(support))
                rm(a_grid_index)
            }
            ## Generate the underlying U grid for the audit
            a_uvec <- sort(c(0, 1, round(rhalton(audit.nu), 8)))
            audit.grid <- list(support = support,
                               uvec = a_uvec,
                               noX = noX)
        } else {
            support <- audit.grid$support
            a_uvec <- audit.grid$uvec
            noX <- audit.grid$noX
        }
      ## Generate the initial constraint grid
        if (noisy) cat("    Generating initial constraint grid...\n")
        if (noX) {
            grid_index <- NULL
        } else {
            initgrid.nx <- min(initgrid.nx, nrow(support))
            grid_index <- sort(
                sample(seq(nrow(support)),
                       initgrid.nx,
                       replace = FALSE,
                       prob = replicate(nrow(support),
                       (1/nrow(support)))))
        }
        if (initgrid.nu > 0) {
            uvec <- sort(c(0, 1, round(rhalton(initgrid.nu), 8)))
        } else {
            uvec <- c(0, 1)
        }
        monoboundAcall <- modcall(call,
                                  newcall = genmonoboundA,
                                  keepargs = monoboundAlist,
                                  newargs = list(m0 = m0,
                                                 m1 = m1,
                                                 uname = uname,
                                                 support = support,
                                                 grid_index = grid_index,
                                                 uvec = uvec,
                                                 splinesobj = splinesobj,
                                                 monov = monov,
                                                 direct = direct))
        ## Generate environment that is to be updated
        modelEnv <- new.env()
        modelEnv$mbobj <- eval(monoboundAcall)
    }
    ## Setup LP problem, or linear components of QCQP problem
    if (!hasArg(equal.coef0) | !hasArg(equal.coef1)) {
        equal.coef0 <- NULL
        equal.coef1 <- NULL
    }
    lpSetup(env = modelEnv, sset = sset, orig.sset = NULL,
            equal.coef0 = equal.coef0, equal.coef1 = equal.coef1,
            solver = solver, direct = direct, rescale = rescale)
    ## Setup QCQP problem
    if (direct) qpSetup(env = modelEnv, sset = sset, rescale = rescale)
    ## Prepare solver messages
    ##
    ## Status codes: 0-unknown; 1-optimal; 2-infeasible; 3-infeasible
    ## or unbounded; 4-unbounded; 5-numerical error; 6-suboptimal;
    ## 7-optimal but infeasible after rescaling; 8-unknown but with a
    ## solution (only for Mosek); 10-maximum iterations reached (only
    ## for Mosek).
    messageAlt <- gsub("\\s+", " ",
                       "If the solver does not return a solution,
                        consider exporting the model
                        and passing it to the solver
                        for more details.\n")
    if (!direct) {
        messageInf <- gsub("\\s+", " ",
                           "Since a minimum criterion was found, the
                        model should be feasible.
                        For more details, consider exporting the model
                        and passing it to the solver.\n")
        messageInfUnb <- gsub("\\s+", " ",
                              "Since a minimum criterion was found, the
                           model is most likely feasible but
                           unbounded. This can happen if
                           the initial grid is too small. Try
                           increasing the parameters 'initgrid.nx'
                           and 'initgrid.nu'. If the model is
                           indeed infeasible,
                           consider exporting the model and passing it
                           to the solver for more details.\n")
    } else {
        messageInf <- gsub("\\s+", " ",
                           "Since the quadratic constraint is constructed
                            without accounting for the shape constraints,
                            infeasibility of the model may arise when the shape
                            constraints are added in the audit procedure. Try
                            increasing the 'criterion.tol' argument to relax the
                            quadratic constraint.\n")
        messageInfUnb <- gsub("\\s+", " ",
                              "The solver may find the model to be infeasible
                               or unbounded if the
                               quadratic constraint is too restrictive
                               (infeasible), or the the initial grid is too
                               small (unbounded).
                               Try increasing the parameters 'criterion.tol'
                               to relax the quadratic constraint; or
                               increasing the parameters 'initgrid.nx'
                               and 'initgrid.nu' to better bound the problem.
                               If the model is
                               indeed infeasible,
                               consider exporting the model and passing it
                               to the solver for more details.\n")
    }
    messageUnb <- gsub("\\s+", " ",
                       "A possible reason for unboundedness is that
                        the initial grid is too small. Try
                        increasing the parameters 'initgrid.nx'
                        and 'initgrid.nu'.\n")
    messageNum <- gsub("\\s+", " ",
                       "A possible reason for numerical issues is that
                        covariates are not scaled appropriately
                        (i.e. the range of magnitudes exceeds 1e13).
                        This is known to cause numerical issues in
                        problems.\n")
    messageSub <- gsub("\\s+", " ",
                       "Tolerance parameters for the solver
                        can be passed through the argument
                        'solver.options'.\n")
    messageOptInf <- gsub("\\s+", " ",
                          "A possible reason for this is that covariates
                           are not scaled appropriately
                           (i.e. the range of magnitudes exceeds 1e13).
                           Tolerance parameters for the solver
                           can also be passed through the argument
                           'solver.options'.\n")
    while (audit_count <= audit.max) {
        if (noisy) {
            cat("\n    Audit count: ", audit_count, "\n", sep = "")
        }
        lpSetupSolver(env = modelEnv, solver = solver)
        if (!direct) {
            lpSetupCriterion(env = modelEnv, sset = sset)
        } else {
            qpSetupCriterion(env = modelEnv)
        }
        minobseq <- criterionMin(env = modelEnv,
                                 sset = sset,
                                 solver = solver,
                                 solver.options = solver.options.criterion,
                                 rescale = rescale,
                                 debug = debug)
        ## Try to diagnose cases where the solution is not
        ## available. This could be due to infeasibility or numerical
        ## issues. To deal with infeasibilty, the LP/QCQP problem is
        ## solved without any shape restrictions. We then check if any
        ## of the lower and upper bounds are violated, which is a
        ## likely cause for infeasible solutions.
        if (minobseq$status %in% c(0, 2, 3, 4, 5)) {
            origMinStatus <- minobseq$status
            ## Stop if issues are numerical, or unbounded, or unknown.
            if (origMinStatus == 0) {
                errMess <-
                    gsub('\\s+', ' ',
                         paste('No solution provided by the solver when
                               minimizing the criterion.',
                               messageAlt))
                stop(errMess, call. = FALSE)
            }
            if (origMinStatus == 4) {
                errMess <-
                    gsub('\\s+', ' ',
                         paste0('No solution to minimizing the criterion
                                since the model is unbounded.',
                                messageUnb))
                stop(errMess, call. = FALSE)
            }
            if (origMinStatus == 5) {
                errMess <-
                    gsub('\\s+', ' ',
                         paste('No solution to minimizing the criterion
                         due to numerical issues.',
                         messageNum))
                stop(errMess, call. = FALSE)
            }
            ## Otherwise, continue and test for infeasibility.
            rm(minobseq)
            if (!direct) {
                lpSetupInfeasible(modelEnv, sset)
            } else {
                qpSetupInfeasible(modelEnv, rescale)
            }
            if (solver == "gurobi") {
                if (debug && solver.options.criterion$outputflag == 1) {
                    cat("Infeasibility diagnosis optimization statistics:\n")
                    cat("------------------------------------------------\n")
                }
            }
            minobseqAlt <- criterionMin(env = modelEnv,
                                        sset = sset,
                                        solver = solver,
                                        solver.options =
                                            solver.options.criterion,
                                        rescale = rescale)
            solVec <- minobseqAlt$x
            ## Test for violation
            negatepos <- which(modelEnv$mbobj$mbs == ">=")
            violateDiff <- modelEnv$mbobj$mbA %*% solVec - modelEnv$mbobj$mbrhs
            violateDiff[negatepos] <- -violateDiff[negatepos]
            violatevec <- as.vector(violateDiff > audit.tol)
            violatepos <- which(violatevec == TRUE)
            violateType <- sapply(violatepos, function(x) {
                if (x %in% modelEnv$mbobj$lb0seq) {
                    if (m0.lb.default == TRUE) {
                        return(paste0("m0.lb = ", round(m0.lb, 6),
                                      " (min. observed outcome by default)"))
                    } else {
                        return(paste0("m0.lb = ", round(m0.lb, 6)))
                    }
                }
                if (x %in% modelEnv$mbobj$lb1seq) {
                    if (m1.lb.default == TRUE) {
                        return(paste0("m1.lb = ", round(m1.lb, 6),
                                      " (min. observed outcome by default)"))
                    } else {
                        return(paste0("m1.lb = ", round(m1.lb, 6)))
                    }
                }
                if (x %in% modelEnv$mbobj$lbteseq) {
                    if (mte.lb.default == TRUE) {
                        return(paste0("mte.lb = ", round(mte.lb, 6),
                                      " (min. treatment effect 6by default)"))
                    } else {
                        return(paste0("mte.lb = ", round(mte.lb, 6)))
                    }
                }
                if (x %in% modelEnv$mbobj$ub0seq) {
                    if (m0.ub.default == TRUE) {
                        return(paste0("m0.ub = ", round(m0.ub, 6),
                                      " (max. observed outcome by default)"))
                    } else {
                        return(paste0("m0.ub = ", round(m0.ub, 6)))
                    }
                }
                if (x %in% modelEnv$mbobj$ub1seq) {
                    if (m1.ub.default == TRUE) {
                        return(paste0("m1.ub = ", round(m1.ub, 6),
                                      " (max. observed outcome by default)"))
                    } else {
                        return(paste0("m1.ub = ", round(m1.ub, 6)))
                    }
                }
                if (x %in% modelEnv$mbobj$ubteseq) {
                    if (mte.ub.default == TRUE) {
                        return(paste0("mte.ub = ", round(mte.ub, 6),
                                      " (max. treatment effect by default)"))
                    } else {
                        return(paste0("mte.ub = ", round(mte.ub, 6)))
                    }
                }
                if (x %in% modelEnv$mbobj$mono0seq) {
                    if (m0.inc == TRUE) return("m0.inc = TRUE")
                    if (m0.dec == TRUE) return("m0.dec = TRUE")
                }
                if (x %in% modelEnv$mbobj$mono1seq) {
                    if (m1.inc == TRUE) return("m1.inc = TRUE")
                    if (m1.dec == TRUE) return("m1.dec = TRUE")
                }
                if (x %in% modelEnv$mbobj$monomteseq) {
                    if (mte.inc == TRUE) return("mte.inc = TRUE")
                    if (mte.dec == TRUE) return("mte.dec = TRUE")
                }
            })
            messageInfDiag <-
                gsub('\\s+', ' ',
                     paste0("The model should only be infeasible if the
                            implied parameter space is empty. The likely cause
                            of an empty parameter space is incoherent shape
                            restrictions. For example, ",
                            paste(unique(violateType), collapse = ", "),
                            " are all set simultaneously. Try changing the
                            shape constraints on the MTR functions.\n"))
            messageUnbDiag <- gsub('\\s+', ' ',
                                   "The model may be unbounded if
                                   the initial grid is too small. Try
                                   increasing the parameters 'initgrid.nx'
                                   and 'initgrid.nu'.\n")
            if (origMinStatus == 2) {
                stop(gsub("\\s+", " ",
                          paste("No solution since the solver proved the
                                 model was infeasible.",
                                messageInfDiag)),
                     call. = FALSE)
            }
            if (origMinStatus == 3) {
                stop(gsub("\\s+", " ",
                          paste("No solution since the solver proved the
                                 model was infeasible or unbounded.",
                                messageInfDiag, messageUnbDiag)),
                     call. = FALSE)
            }
        }
        if (!direct && noisy) {
            cat("    Minimum criterion: ", fmtResult(minobseq$obj), "\n",
                sep = "")
        }
        if (direct && noisy) {
            cat("    Minimum criterion: ",
                fmtResult(minobseq$obj * modelEnv$drN + modelEnv$ssy),
                "\n",
                sep = "")
        }
        ## Perform specification test
        if (!direct && !is.null(orig.sset) && !is.null(orig.criterion)) {
            lpSetupCriterionBoot(modelEnv, sset, orig.sset,
                                 orig.criterion, criterion.tol,
                                 setup = TRUE)
            minobseqTest <- criterionMin(env = modelEnv,
                                         sset = sset,
                                         solver = solver,
                                         solver.options =
                                             solver.options.criterion,
                                         rescale = rescale)
            lpSetupCriterionBoot(modelEnv, sset, orig.sset,
                                 orig.criterion, criterion.tol,
                                 setup = FALSE)
        }
        ## Provide warnings if solutions are suboptimal.
        bWarn <- NULL
        if (minobseq$status == 6) {
            bWarn <-
                paste(bWarn,
                      gsub("\\s+", " ",
                           paste('The solver was unable to satisfy
                               the optimality tolerance when
                               minimizing the criterion, so a suboptimal
                               solution is returned.')))
            bWarn <- paste(bWarn, messageSub)
        }
        if (minobseq$status == 7) {
            bWarn <-
                paste(bWarn,
                      gsub("\\s+", " ",
                           paste('The solution to the problem of
                                  minimizing the criterion
                                  is optimal, but infeasible after
                                  rescaling.')))
            bWarn <- paste(bWarn, messageOptInf)
        }
        if (minobseq$status == 8) {
            bWarn <-
                paste(bWarn,
                      gsub("\\s+", " ",
                           paste("The solution status (e.g. 'OPTIMAL')
                                  to the problem of
                                  minimizing the criterion
                                  is unknown---Rmosek provided a solution
                                  but did not provide a status.")))
            bWarn <- paste(bWarn, messageOptInf)
        }
        if (minobseq$status == 10) {
            bWarn <-
                paste(bWarn,
                      gsub("\\s+", " ",
                           paste("Minimizing the criterion was terminated
                                  by Rmosek because the maximum iteration limit
                                  was reached.")))
            bWarn <- paste(bWarn, messageOptInf)
        }
        if (!is.null(bWarn)) warning(bWarn, call. = FALSE,
                                     immediate. = TRUE)

        ## Obtain bounds
        if (noisy) {
            cat("    Obtaining bounds...\n")
        }
        if (!direct) {
            lpSetupBound(env = modelEnv,
                         g0 = gstar0,
                         g1 = gstar1,
                         sset = sset,
                         criterion.tol = criterion.tol,
                         criterion.min = minobseq$obj,
                         solver = solver,
                         setup = TRUE)
        } else {
            qpSetupBound(env = modelEnv,
                         g0 = gstar0,
                         g1 = gstar1,
                         criterion.tol = criterion.tol,
                         criterion.min = minobseq$obj,
                         rescale = rescale)
        }
        result <- bound(env = modelEnv,
                          sset = sset,
                          solver = solver,
                          solver.options = solver.options.bounds,
                          noisy = noisy,
                          smallreturnlist = smallreturnlist,
                          rescale = rescale,
                          debug = debug)
        if (result$error == TRUE) {
            errMess <- NULL
            errTypes <- NULL
            for (type in c('min', 'max')) {
                tmpName <- paste0(type, 'status')
                if (type == 'min') tmpType <- 'minimization'
                if (type == 'max') tmpType <- 'maximization'
                if (result[[tmpName]] == 0) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The solver did not provide a
                                     solution for the', tmpType,
                                     'problem.')))
                }
                if (result[[tmpName]] == 2) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType, 'problem was proven
                                          to be infeasible.')))
                }
                if (result[[tmpName]] == 3) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType,
                                         'problem was proven to be
                                          infeasible or unbounded.')))
                }
                if (result[[tmpName]] == 4) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType,
                                         'problem was proven to be
                                          unbounded.')))
                }
                if (result[[tmpName]] == 5) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType,
                                         'problem resulted in
                                          numerical issues.')))
                }
                errTypes <- c(errTypes, result[[tmpName]])
            }
            ## Include explanation
            origErrTypes <- sort(unique(errTypes))
            errTypes <- sort(unique(errTypes))
            if (length(errTypes) == 2 &&
                errTypes[1] == 2 && errTypes[2] == 3) errTypes <- 3
            for (et in errTypes) {
                if (et == 0) {
                    errMess <- paste(errMess, messageAlt)
                }
                if (et == 2) {
                    errMess <- paste(errMess, messageInf)
                }
                if (et == 3) {
                    errMess <- paste(errMess, messageInfUnb)
                }
                if (et == 4) {
                    errMess <- paste(errMess, messageUnb)
                }
                if (et == 5) {
                    errMess <- paste(errMess, messageNum)
                }
            }
            return(list(error = errMess,
                        errorTypes = origErrTypes,
                        min = result$min,
                        max = result$max,
                        status.min = result$minstatus,
                        status.max = result$maxstatus,
                        audit.criterion = minobseq$obj * modelEnv$drN +
                            modelEnv$ssy,
                        audit.criterion.raw = minobseq$obj,
                        audit.criterion.status = minobseq$status,
                        audit.count = audit_count - 1,
                        audit.grid = audit.grid))
        }
        ## Provide warnings if solutions are suboptimal.
        bWarn <- NULL
        bWarnTypes <- NULL
        for (type in c('min', 'max')) {
            tmpName <- paste0(type, 'status')
            if (type == 'min') tmpType <- 'minimization'
            if (type == 'max') tmpType <- 'maximization'
            if (result[[tmpName]] == 6) {
                bWarn <-
                    paste(bWarn,
                          gsub("\\s+", " ",
                               paste('The solver was unable to satisfy
                               the optimality tolerance for the',
                               tmpType, 'problem, so a suboptimal
                               solution is returned.')))
            }
            if (result[[tmpName]] == 7) {
                bWarn <-
                    paste(bWarn,
                          gsub("\\s+", " ",
                               paste('The solution to the',
                                     tmpType, 'problem is optimal,
                                     but infeasible after rescaling.')))
            }
            if (result[[tmpName]] == 8) {
                bWarn <-
                    paste(bWarn,
                          gsub("\\s+", " ",
                               paste("Rmosek did not provide a solution
                                          status (e.g. 'OPTIMAL') to the",
                                     tmpType, 'problem.')))
            }
            if (result[[tmpName]] == 10) {
                bWarn <-
                    paste(bWarn,
                          gsub("\\s+", " ",
                               paste("Rmosek terminated the",
                                     tmpType, 'problem because the
                                         iteration limit was reached.')))
            }
            bWarnTypes <- c(bWarnTypes,
                            result[[tmpName]])
        }
        bWarnTypes <- sort(unique(bWarnTypes))
        for (wt in bWarnTypes) {
            if (wt == 6) {
                bWarn <- paste(bWarn, messageSub)
            }
            if (wt == 7) {
                bWarn <- paste(bWarn, messageOptInf)
            }
        }
        if (!is.null(bWarn)) warning(bWarn, call. = FALSE, immediate. = TRUE)
        ## Save results
        solVecMin <- c(result$ming0, result$ming1)
        solVecMax <- c(result$maxg0, result$maxg1)
        optstatus <- min(c(result$minstatus,
                           result$maxstatus))
        if (existsolution == FALSE) existsolution <- TRUE
        prevbound <- c(result$min, result$max)

        ## Test for violations of shape constraints when obtaining the bounds
        monoboundAcall <- modcall(call,
                                  newcall = genmonoboundA,
                                  keepargs = monoboundAlist,
                                  newargs = list(m0 = m0,
                                                 m1 = m1,
                                                 uname = uname,
                                                 support = support,
                                                 grid_index =
                                                     seq(nrow(support)),
                                                 uvec = a_uvec,
                                                 splinesobj =
                                                     splinesobj,
                                                 monov = monov,
                                                 solution.m0.min =
                                                     result$ming0,
                                                 solution.m1.min =
                                                     result$ming1,
                                                 solution.m0.max =
                                                     result$maxg0,
                                                 solution.m1.max =
                                                     result$maxg1,
                                                 audit.tol = audit.tol,
                                                 direct = direct))
        auditObj <- eval(monoboundAcall)
        ## Combine the violation matrices
        violateMat <- NULL
        if (!is.null(auditObj$bounds)) {
            violateMat <- rbind(violateMat, auditObj$bounds$violateMat)
            auditObj$bounds$violateMat <- NULL
        }
        if (!is.null(auditObj$mono)) {
            violateMat <- rbind(violateMat, auditObj$mono$violateMat)
            auditObj$mono$violateMat <- NULL
        }
        ## Deal with possible violations when audit and initial grid match
        if (initgrid.nx == audit.nx &&
            initgrid.nu == audit.nu) {
            if (!is.null(violateMat) && nrow(violateMat) > 0) {
                origTol <- audit.tol
                origViolations <- nrow(violateMat)
                origViolationsMagnitude <- max(violateMat$diff)
                ## If violations occur when initial grid and audit grid
                ## match, then there is a precision issue. By default, the
                ## tolerance in the package will be adjusted to be 'twice'
                ## the magnitude of solver's tolerance.
                audit.tol <- (audit.tol / (10^magnitude(audit.tol))) *
                    (10^(magnitude(audit.tol) / 2))
                violateMat <- violateMat[violateMat$diff > audit.tol, ]
            } else {
                origViolations <- 0
            }
            ## If the audit grid and initial grid are identical, then
            ## the audit must end after the first iteration---there
            ## are no points to add.
            if (!is.null(violateMat) && nrow(violateMat) > 0) {
                warning(paste0(gsub("\\s+", " ",
                                    paste0("Audit finished: violations
                                           continue to occur although
                                           the constraint grid and
                                           audit grid are identical,
                                           and the audit tolerance has
                                           increased from ",
                                           format(origTol, scientific = TRUE),
                                           " to ",
                                           format(audit.tol, scientific = TRUE),
                                           ". This suggests precision of the
                                           solver exceeds that of R,
                                           and may be due to scaling issues
                                           in the model.
                                           Setting 'audit.tol' to be larger
                                           than the maximum violation of ",
                                           format(origViolationsMagnitude,
                                                  scientific = TRUE),
                                           " will eliminate this message,
                                           but will not address the nature of
                                           the problem.")),
                               "\n"),
                        call. = FALSE, immediate. = TRUE)
                break
            } else {
                if (origViolations > 0) {
                    warning(gsub("\\s+", " ",
                                 paste0(origViolations,
                                        " violations originally found despite
                                         the constraint grid and audit grid
                                         being identical.
                                    Audit tolerance was increased from ",
                                    format(origTol, scientific = TRUE), " to ",
                                    format(audit.tol, scientific = TRUE),
                                    ".  This suggests precision of the
                                 solver exceeds that of R,
                                 and may be due to scaling issues
                                 in the model.")), "\n",
                            call. = FALSE,
                            immediate. = TRUE)
                }
                break
            }
        }
        ## Address violations by expanding initial grid
        if (!is.null(violateMat)) {
            if (noisy) cat("    Violations: ", nrow(violateMat), "\n")
            if (audit_count == audit.max) {
                ## End audit and return violation matrix if max audit
                ## is achieved.
                if (noisy) {
                    warning(paste0(gsub("\\s+", " ",
                                        paste0("Audit finished: maximum number
                                                of audits (audit.max = ",
                                               audit.max,
                                               ") reached. Try increasing
                                                audit.max.")), "\n"),
                            call. = FALSE, immediate. = TRUE)
                }
                break
            } else {
                ## Expand the initial grid if audits are still possible.
                ## If number of violations is fewer than audit.add, simply
                ## include all points
                if (nrow(violateMat) > audit.add) {
                    ## If number of violations exceeds audit.add, then a
                    ## selection rule is applied. Sort the violations by
                    ## type and size.
                    violateMat <- violateMat[order(violateMat$type,
                                                   violateMat$grid.x,
                                                   -violateMat$diff), ]
                    violateMat$group.count <-
                        unlist(lapply(table(violateMat$group.name),
                                      function(x) seq(x)))
                    if (nrow(violateMat[violateMat$group.count <= 1, ]) >=
                        audit.add) {
                        violateMat <- violateMat[violateMat$group.count <= 1, ]
                    } else {
                        k <- 2
                        for (i in 2:max(violateMat$group.count)) {
                            if (nrow(violateMat[violateMat$group.count <= k, ])
                                >= audit.add) {
                                break
                            } else {
                                k <- k + 1
                            }
                        }
                        vmatTmp <- violateMat[violateMat$group.count <=
                                              (k - 1), ]
                        vmatAdd <- violateMat[violateMat$group.count == k, ]
                        vmatAdd <- vmatAdd[order(-vmatAdd$diff), ]
                        vmatAdd <- vmatAdd[1:(audit.add - nrow(vmatTmp)), ]
                        violateMat <- rbind(vmatTmp, vmatAdd)
                        violateMat <- violateMat[order(violateMat$type,
                                                       violateMat$grid.x,
                                                       -violateMat$diff), ]
                    }
                }
                ## Now expand the initial grid
                if (noisy) {
                    if (nrow(violateMat) > 1) ps <- 'points'
                    if (nrow(violateMat) == 1) ps <- 'point'
                    cat("    ",
                        gsub("\\s+", " ",
                             paste0("Expanding constraint grid to
                                        include ", nrow(violateMat),
                                    " additional ", ps, "...")), "\n", sep = "")
                }
                types <- unique(violateMat$type)
                addm0 <- NULL
                addm1 <- NULL
                addmte <- NULL
                ## Placeholders for violation types
                tmpAdd <- 0
                addlb0seq <- NULL
                addlb1seq <- NULL
                addlbteseq <- NULL
                addub0seq <- NULL
                addub1seq <- NULL
                addubteseq <- NULL
                addmono0incseq <- NULL
                addmono0decseq <- NULL
                addmono1incseq <- NULL
                addmono1decseq <- NULL
                addmonoteincseq <- NULL
                addmonotedecseq <- NULL
                nShapeConstraints <- max(c(modelEnv$mbobj$lb0seq,
                                           modelEnv$mbobj$lb1seq,
                                           modelEnv$mbobj$ub0seq,
                                           modelEnv$mbobj$ub1seq,
                                           modelEnv$mbobj$lbteseq,
                                           modelEnv$mbobj$ubteseq,
                                           modelEnv$mbobj$mono0seq[, 1],
                                           modelEnv$mbobj$mono1seq[, 1],
                                           modelEnv$mbobj$monoteseq[, 1]))
                for (i in types) {
                    ## Expand constraints for m0
                    if (i == 1) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        addlb0seq <- seq(length(addIndex)) + tmpAdd
                        tmpAdd <- tmpAdd + length(addIndex)
                        tmpGrid <- auditObj$bounds$bdA$m0.lb
                        modelEnv$model$sense <- c(modelEnv$model$sense,
                                               rep(">=", length(addIndex)))
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(m0.lb,
                                                 length(addIndex)))
                        addm0 <- rbind(addm0, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        auditObj$bounds$bdA$m0.lb <- NULL
                    }
                    if (i == 4) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        addub0seq <- seq(length(addIndex)) + tmpAdd
                        tmpAdd <- tmpAdd + length(addIndex)
                        tmpGrid <- auditObj$bounds$bdA$m0.ub
                        modelEnv$model$sense <- c(modelEnv$model$sense,
                                               rep("<=", length(addIndex)))
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(m0.ub,
                                                 length(addIndex)))
                        addm0 <- rbind(addm0, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        auditObj$bounds$bdA$m0.ub <- NULL
                    }
                    if (i %in% c(7, 8)) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        tmpGrid <- auditObj$mono$monoA$m0
                        if (i == 7) {
                            addmono0incseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            modelEnv$model$sense <- c(modelEnv$model$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
                            addmono0decseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            modelEnv$model$sense <- c(modelEnv$model$sense,
                                                   rep("<=", length(addIndex)))
                        }
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(0, length(addIndex)))
                        addm0 <- rbind(addm0, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        if (i == 8) auditObj$mono$monoA$m0 <- NULL
                    }
                }
                if (!is.null(addm0)) {
                    ## Update the constraint matrix
                    if (!direct) addCol <- length(sset$s1$g1)
                    if (direct) addCol <- ncol(sset$s1$g1)
                    tmpMat <- cbind(matrix(0, nrow = nrow(addm0),
                                           ncol = 2 * sn),
                                    addm0,
                                    matrix(0, nrow = nrow(addm0),
                                           ncol = addCol))
                    if (rescale) {
                        tmpMat <- sweep(x = tmpMat,
                                        MARGIN = 2,
                                        STATS = modelEnv$colNorms,
                                        FUN = '/')
                    }
                    modelEnv$model$A <- rbind(modelEnv$model$A, tmpMat)
                    rm(addCol, tmpMat)
                    ## Update the contraint sequences
                    modelEnv$mbobj$lb0seq <- c(modelEnv$mbobj$lb0seq,
                                            addlb0seq + nShapeConstraints)
                    rm(addlb0seq)
                    modelEnv$mbobj$ub0seq <- c(modelEnv$mbobj$ub0seq,
                                            addub0seq + nShapeConstraints)
                    rm(addub0seq)
                    if (!is.null(addmono0incseq)) {
                        modelEnv$mbobj$mono0seq <-
                            rbind(modelEnv$mbobj$mono0seq,
                                  cbind(addmono0incseq + nShapeConstraints, 1))
                    }
                    if (!is.null(addmono0decseq)) {
                        modelEnv$mbobj$mono0seq <-
                            rbind(modelEnv$mbobj$mono0seq,
                                  cbind(addmono0decseq + nShapeConstraints, -1))
                    }
                    rm(addmono0incseq, addmono0decseq)
                }
                rm(addm0)
                ## Expand constraints for m1
                for (i in types) {
                    if (i == 2) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        addlb1seq <- seq(length(addIndex)) + tmpAdd
                        tmpAdd <- tmpAdd + length(addIndex)
                        tmpGrid <- auditObj$bounds$bdA$m1.lb
                        modelEnv$model$sense <- c(modelEnv$model$sense,
                                               rep(">=", length(addIndex)))
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(m1.lb,
                                                 length(addIndex)))
                        addm1 <- rbind(addm1, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        auditObj$bounds$bdA$m1.lb <- NULL
                    }
                    if (i == 5) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        addub1seq <- seq(length(addIndex)) + tmpAdd
                        tmpAdd <- tmpAdd + length(addIndex)
                        tmpGrid <- auditObj$bounds$bdA$m1.ub
                        modelEnv$model$sense <- c(modelEnv$model$sense,
                                               rep("<=", length(addIndex)))
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(m1.ub,
                                                 length(addIndex)))
                        addm1 <- rbind(addm1, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        auditObj$bounds$bdA$m1.ub <- NULL
                    }
                    if (i %in% c(9, 10)) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        tmpGrid <- auditObj$mono$monoA$m1
                        if (i == 9) {
                            addmono1incseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            modelEnv$model$sense <- c(modelEnv$model$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
                            addmono1decseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            modelEnv$model$sense <- c(modelEnv$model$sense,
                                                   rep("<=", length(addIndex)))
                        }
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(0, length(addIndex)))
                        addm1 <- rbind(addm1, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        if (i == 10) auditObj$mono$monoA$m1 <- NULL
                    }
                }
                if (!is.null(addm1)) {
                    ## Update the constraint matrix
                    if (!direct) addCol <- length(sset$s1$g0)
                    if (direct) addCol <- ncol(sset$s1$g0)
                    tmpMat <- cbind(matrix(0, nrow = nrow(addm1),
                                           ncol = 2 * sn),
                                    matrix(0, nrow = nrow(addm1),
                                           ncol = addCol),
                                    addm1)
                    if (rescale) {
                        tmpMat <- sweep(x = tmpMat,
                                        MARGIN = 2,
                                        STATS = modelEnv$colNorms,
                                        FUN = '/')
                    }
                    modelEnv$model$A <-
                        rbind(modelEnv$model$A, tmpMat)
                    rm(addCol, tmpMat)
                    ## Update the contraint sequences
                    modelEnv$mbobj$lb1seq <- c(modelEnv$mbobj$lb1seq,
                                            addlb1seq + nShapeConstraints)
                    rm(addlb1seq)
                    modelEnv$mbobj$ub1seq <- c(modelEnv$mbobj$ub1seq,
                                            addub1seq + nShapeConstraints)
                    rm(addub1seq)
                    if (!is.null(addmono1incseq)) {
                        modelEnv$mbobj$mono1seq <-
                            rbind(modelEnv$mbobj$mono1seq,
                                  cbind(addmono1incseq + nShapeConstraints, 1))
                    }
                    if (!is.null(addmono1decseq)) {
                        modelEnv$mbobj$mono1seq <-
                            rbind(modelEnv$mbobj$mono1seq,
                                  cbind(addmono1decseq + nShapeConstraints, -1))
                    }
                    rm(addmono1incseq, addmono1decseq)
                }
                rm(addm1)
                ## Expand constraints for mte
                for (i in types) {
                    if (i == 3) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        addlbteseq <- seq(length(addIndex)) + tmpAdd
                        tmpAdd <- tmpAdd + length(addIndex)
                        tmpGrid <- auditObj$bounds$bdA$mte.lb
                        modelEnv$model$sense <- c(modelEnv$model$sense,
                                               rep(">=", length(addIndex)))
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(mte.lb,
                                                 length(addIndex)))
                        addmte <- rbind(addmte, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        auditObj$bounds$bdA$mte.lb <- NULL
                    }
                    if (i == 6) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        addubteseq <- seq(length(addIndex)) + tmpAdd
                        tmpAdd <- tmpAdd + length(addIndex)
                        tmpGrid <- auditObj$bounds$bdA$mte.ub
                        modelEnv$model$sense <- c(modelEnv$model$sense,
                                               rep("<=", length(addIndex)))
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(mte.ub,
                                                 length(addIndex)))
                        addmte <- rbind(addmte, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        auditObj$bounds$bdA$mte.lb <- NULL
                    }
                    if (i %in% c(11, 12)) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        tmpGrid <- auditObj$mono$monoA$mte
                        if (i == 11) {
                            addmonoteincseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            modelEnv$model$sense <- c(modelEnv$model$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
                            addmonotedecseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            modelEnv$model$sense <- c(modelEnv$model$sense,
                                                   rep("<=", length(addIndex)))
                        }
                        modelEnv$model$rhs <- c(modelEnv$model$rhs,
                                             rep(0, length(addIndex)))
                        addmte <- rbind(addmte, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        if (i == 12) auditObj$mono$monoA$mte <- NULL
                    }
                }
                if (!is.null(addmte)) {
                    ## Update the constraint matrix
                    tmpMat <- cbind(matrix(0, nrow = nrow(addmte),
                                           ncol = 2 * sn),
                                    addmte)
                    if (rescale) {
                        tmpMat <- sweep(x = tmpMat,
                                        MARGIN = 2,
                                        STATS = modelEnv$colNorms,
                                        FUN = '/')
                    }
                    modelEnv$model$A <-
                        rbind(modelEnv$model$A, tmpMat)
                    rm(tmpMat)
                    ## Update the contraint sequences
                    modelEnv$mbobj$lbteseq <- c(modelEnv$mbobj$lbteseq,
                                             addlbteseq + nShapeConstraints)
                    rm(addlbteseq)
                    modelEnv$mbobj$ubteseq <- c(modelEnv$mbobj$ubteseq,
                                             addubteseq + nShapeConstraints)
                    rm(addubteseq)
                    if (!is.null(addmonoteincseq)) {
                        modelEnv$mbobj$monoteseq <-
                            rbind(modelEnv$mbobj$monoteseq,
                                  cbind(addmonoteincseq + nShapeConstraints, 1))
                    }
                    if (!is.null(addmonotedecseq)) {
                        modelEnv$mbobj$monoteseq <-
                            rbind(modelEnv$mbobj$monoteseq,
                                  cbind(addmonotedecseq + nShapeConstraints,
                                        -1))
                    }
                    rm(addmonoteincseq, addmonotedecseq)
                }
                rm(addmte)
                ## Move on to next iteration of the audit
                audit_count <- audit_count + 1
                if (!direct) lpSetupBound(env = modelEnv, setup = FALSE)
                if (direct)  qpSetupBound(env = modelEnv, setup = FALSE)
            }
        } else {
            ## If no violations, then end the audit
            if (noisy) {
                cat("    Violations: 0\n")
                cat("    Audit finished.\n\n")
            }
            break
        }
    }
    if (!is.null(violateMat)) {
        violateMat$type <- violateMat$type.string
        violateMat$type.string <- NULL
        violateMat$group.name <- NULL
        violateMat$pos <- NULL
    }
    ## Clean up status codes
    result[['minstatus']] <- statusString(result[['minstatus']], solver)
    result[['maxstatus']] <- statusString(result[['maxstatus']], solver)
    minobseq$status <- statusString(minobseq$status, solver)
    ## Return output
    output <- list(max = result$max,
                   min = result$min,
                   result = result,
                   gridobj = list(audit.grid = audit.grid,
                                  violations = violateMat),
                   auditcount = audit_count,
                   minobseq = minobseq$obj,
                   minobseq.status = minobseq$status)
    if (direct) {
        output$minobseq <- output$minobseq * modelEnv$drN + modelEnv$ssy
        output$minobseq.raw <- minobseq$obj
    }
    if (!is.null(orig.sset) && !is.null(orig.criterion)) {
        output$spectest = minobseqTest$obj
    }
    return(output)
}

#' Select points from audit grid to add to the constraint grid
#'
#' This function selects which points from the audit grid should be
#' included into the original grid. Both the constraint grid and audit
#' grid are represented as constraints in an LP/QCQP problem. This
#' function selects which points in the audit grid (i.e. which rows in
#' the audit constraint matrix) should be added to the constraint grid
#' (i.e. should be appended to the constraint matrix).
#'
#' @param diffVec numeric vector, with a positive value indicating a
#'     violation of a shape constraint.
#' @param audit.add integer, the number of points from the audit grid
#'     to add to the initial for each constraint type. For instance, if
#'     there are 5 different kinds of constraints imposed, and
#'     \code{audit.add = 5}, then up to 30 points may be added to the
#'     constraint grid.
#' @param lb0seq integer vector, indicates which rows in the audit
#'     constraint matrix correspond to the lower bound for m0.
#' @param lb1seq integer vector, indicates which rows in the audit
#'     constraint matrix correspond to the lower bound for m1.
#' @param lbteseq integer vector, indicates which rows in the audit
#'     constriant matrix correspond to the lower bound for the
#'     treatment effect.
#' @param ub0seq integer vector, indicates which rows in the audit
#'     constraint matrix correspond to the upper bound for m0.
#' @param ub1seq integer vector, indicates which rows in the audit
#'     constraint matrix correspond to the upper bound for m1.
#' @param ubteseq integer vector, indicates which rows in the audit
#'     constriant matrix correspond to the upper bound for the
#'     treatment effect.
#' @param mono0seq integer matrix, indicates which rows in the audit
#'     constraint matrix correspond to the monotonicity conditions for
#'     m0, and whether the constraint is increasing (+1) or decreasing
#'     (-1).
#' @param mono1seq integer matrix, indicates which rows in the audit
#'     constraint matrix correspond to the monotonicity conditions for
#'     m1, and whether the constraint is increasing (+1) or decreasing
#'     (-1).
#' @param monoteseq integer matrix, indicates which rows in the audit
#'     constraint matrix correspond to the monotonicity conditions for
#'     the treatment effect, and whether the constraint is increasing
#'     (+1) or decreasing (-1).
#' @param mbmap integer vector, indexes the X-value associated with
#'     each row in the audit constraint matrix.
#' @return The audit grid is represented using a set of constraint
#'     matrices. Each point in the audit grid corresponds to a set of
#'     rows in the constraint matrices. The function simply returns
#'     the vector of row numbers for the points from the audit grid
#'     whose corresponding constraints should be added to the original
#'     LP/QCQP problem (i.e. the points to add to the original grid).
selectViolations <- function(diffVec, audit.add,
                             lb0seq, lb1seq, lbteseq,
                             ub0seq, ub1seq, ubteseq,
                             mono0seq, mono1seq, monoteseq,
                             mbmap) {
    typeVec <- c(rep(1, times = length(lb0seq)),
                 rep(2, times = length(lb1seq)),
                 rep(3, times = length(lbteseq)),
                 rep(4, times = length(ub0seq)),
                 rep(5, times = length(ub1seq)),
                 rep(6, times = length(ubteseq)),
                 rep(7, sum(times = mono0seq[, 2] > 0)),
                 rep(8, sum(times = mono0seq[, 2] < 0)),
                 rep(9, sum(times = mono1seq[, 2] > 0)),
                 rep(10, sum(times = mono1seq[, 2] < 0)),
                 rep(11, sum(times = monoteseq[, 2] > 0)),
                 rep(12, sum(times = monoteseq[, 2] < 0)))
    ## Store all points that violate the constraints
    violateMat <- data.frame(cbind(c(lb0seq, lb1seq, lbteseq,
                                     ub0seq, ub1seq, ubteseq,
                                     mono0seq[, 1],
                                     mono1seq[, 1],
                                     monoteseq[, 1]),
                                   typeVec,
                                   mbmap))
    colnames(violateMat) <- c("row", "type", "grid.x")
    violateMat$diff <- diffVec
    violateMat <- violateMat[violateMat$diff > 0, ]
    if (nrow(violateMat) < audit.add) {
        return(violateMat)
    } else {
        ## For each point in the X-grid, find the U that violats
        ## the constraints the most
        violateMat <- violateMat[order(violateMat$type,
                                       violateMat$grid.x,
                                       -violateMat$diff), ]
        violateMat$group.name <- paste0(violateMat$type,
                                        ".", violateMat$grid.x)
        violateMat$group.count <- unlist(lapply(table(violateMat$group.name),
                                                function(x) seq(x)))
        if (nrow(violateMat[violateMat$group.count <= 1, ]) >= audit.add) {
            return(violateMat[violateMat$group.count <= 1, ])
        } else {
            k <- 2
            for (i in 2:max(violateMat$group.count)) {
                if (nrow(violateMat[violateMat$group.count <= k, ]) >=
                    audit.add) {
                    break
                } else {
                    k <- k + 1
                }
            }
            vmatTmp <- violateMat[violateMat$group.count <= (k - 1), ]
            vmatAdd <- violateMat[violateMat$group.count == k, ]
            vmatAdd <- vmatAdd[order(-vmatAdd$diff), ]
            vmatAdd <- vmatAdd[1:(audit.add - nrow(vmatTmp)), ]
            violateMat <- rbind(vmatTmp, vmatAdd)
            violateMat <- violateMat[order(violateMat$type,
                                           violateMat$grid.x,
                                           -violateMat$diff), ]
            return(violateMat)
        }
    }
}

#' Generate Halton sequence
#'
#' This function generates a one dimensional Halton sequence.
#'
#' @param n Number of draws.
#' @param base Base used for the Halton sequence, set to 2 by default.
#' @return A sequence of randomly drawn numbers.
#' @export
rhalton <- function(n, base = 2) {
    output <- NULL
    for (j in 1:n) {
        f <- 1
        r <- 0
        i <- j
        while(i > 0) {
            f  <- f / base
            r <- r + f * (i %% base)
            i <- floor(i / base)
        }
        output <- c(output, r)
    }
    return(output)
}

#' Convert status code to string
#'
#' This function returns the status code specific to a solver.
#'
#' @param status Status code.
#' @param solver Name of solver, either 'gurobi', 'cplexapi', or
#'     'lpsolveapi'.
#' @return Status specific to solver, e.g. 'OPTIMAL (2)'.
statusString <- function(status, solver) {
    if (solver == 'gurobi') {
        if (status == 1) statusStr <- 'OPTIMAL (2)'
        if (status == 2) statusStr <- 'INFEASIBLE (3)'
        if (status == 3) statusStr <- 'INF_OR_UNBD (4)'
        if (status == 4) statusStr <- 'UNBOUNDED (5)'
        if (status == 5) statusStr <- 'NUMERIC (12)'
        if (status == 6) statusStr <- 'SUBOPTIMAL (13)'
    }
    if (solver == 'cplexapi') {
        if (status == 1) statusStr <-
                             'CPX_STAT_OPTIMAL (1)'
        if (status == 2) statusStr <-
                             'CPX_STAT_INFEASIBLE (3)'
        if (status == 3) statusStr <-
                             'CPX_STAT_INForUNBD (4)'
        if (status == 4) statusStr <-
                             'CPX_STAT_UNBOUNDED (2)'
        if (status == 6) statusStr <-
                             'CPX_STAT_NUM_BEST (6)'
        if (status == 7) statusStr <-
                             'CPX_STAT_OPTIMAL_INFEAS (5)'
    }
    if (solver == 'lpsolveapi') {
        if (status == 1) statusStr <- 'Optimal (0)'
        if (status == 2) statusStr <- 'Infeasible (2)'
        if (status == 4) statusStr <- 'Unbounded (3)'
        if (status == 5) statusStr <-
                             'Numerical failure (5)'
        if (status == 6) statusStr <- 'Sub-optimal (1)'
    }
    if (solver == 'rmosek') {
        if (status == 1) statusStr <- 'OPTIMAL'
        if (status == 4) statusStr <- 'DUAL_INFEASIBLE_CER'
        if (status == 2) statusStr <- 'PRIMAL_INFEASIBLE_CER'
        if (status == 5) statusStr <- 'NUMERICAL_PROBLEM (10025)'
        if (status == 6) statusStr <- 'STALL (10006)'
        if (status == 8) statusStr <- 'UNKNOWN'
        if (status == 10) statusStr <- 'MAX_ITERATIONS (10000)'
    }
    if (status == 0) statusStr <- 'Unknown error'
    return(statusStr)
}
