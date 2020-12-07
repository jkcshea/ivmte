#' Audit procedure
#'
#' This is the wrapper for running the entire audit procedure. This
#' function sets up the LP problem of minimizing the violation of
#' observational equivalence for the set of IV-like estimands, while
#' satisfying boundedness and monotonicity constraints declared by the
#' user. Rather than enforce that boundedness and monotonicity hold
#' across the entire support of covariates and unobservables, this
#' procedure enforces the conditions over a grid of points. This grid
#' corresponds to the set of values the covariates can take, and a set
#' of values of the unobservable term. The size of this grid is
#' specified by the user in the function arguments. The procedure
#' first estimates the bounds while imposing the shape constraints for
#' an initial subset of points in the grid. The procedure then goes on
#' to check ('audit') whether the constraints are satisfied over the
#' entire grid. Any point where either the boundedness or monotonicity
#' constraints are violated are incorporated into the initial grid,
#' and the process is repeated until the audit no longer finds any
#' violations, or until some maximum number of iterations is
#' reached.
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
#'
#' @inheritParams ivmteEstimate
#'
#' @return a list. Included in the list are estimates of the treatment
#'     effect bounds; the minimum violation of observational
#'     equivalence of the set of IV-like estimands; the list of
#'     matrices and vectors defining the LP problem; the points used
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
#'       lpsolver = "lpSolveAPI")
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
                  sset, gstar0, gstar1,
                  orig.sset = NULL, orig.criterion = NULL,
                  criterion.tol = 0,
                  lpsolver, lpsolver.options, lpsolver.presolve,
                  lpsolver.options.criterion, lpsolver.options.bounds,
                  smallreturnlist = FALSE,
                  noisy = TRUE, debug = FALSE) {
    call  <- match.call()
    lpsolver <- tolower(lpsolver)
    ## Determine if whether IV-like moments or direct MTR regression
    ## will be used.
    direct <- ('Q' %in% names(sset$s1))
    ## Set the audit tolerance
    if (!hasArg(audit.tol)) {
        if (lpsolver == "gurobi") {
            if (hasArg(lpsolver.options) &&
                !is.null(lpsolver.options$FeasibilityTol)) {
                audit.tol <- lpsolver.options$FeasibilityTol
            } else if (hasArg(lpsolver.options.bounds) &&
                       !is.null(lpsolver.options.bounds$FeasibilityTol)) {
                audit.tol <- lpsolver.options.bounds$FeasibilityTol
            } else {
                audit.tol <- 1e-06
            }
        }
        if (lpsolver == "cplexapi") {
            if (hasArg(lpsolver.options)) {
                audit.tol <- optionsCplexAPITol(lpsolver.options)
            } else if (hasArg(lpsolver.options.bounds)){
                audit.tol <- optionsCplexAPITol(lpsolver.options.bounds)
            } else {
                audit.tol <- 1e-06
            }
        }
        if (lpsolver == "lpsolveapi") {
            audit.tol <- 1e-06
        }
    }
    ## Organize LP options
    if (lpsolver == "gurobi") {
        ## Construct default options
        lpsolver.options.default <- list(dualreductions = 1,
                                         FeasibilityTol = 1e-06)
        if (debug)  lpsolver.options.default$outputflag <- 1
        if (!debug) lpsolver.options.default$outputflag <- 0
        if (hasArg(lpsolver.presolve)) lpsolver.options.default$presolve <-
                                           as.integer(lpsolver.presolve)
        ## Prepare user options
        if (hasArg(lpsolver.options)) {
            lpsolver.options <- optionsGurobi(lpsolver.options, debug)
            if (hasArg(lpsolver.presolve)) lpsolver.options$presolve <-
                                               as.integer(lpsolver.presolve)
            lpsolver.options.criterion <- lpsolver.options
            lpsolver.options.bounds <- lpsolver.options
        } else {
            if (hasArg(lpsolver.options.criterion)) {
                lpsolver.options.criterion <-
                    optionsGurobi(lpsolver.options.criterion, debug)
                if (hasArg(lpsolver.presolve)) {
                    lpsolver.options.criterion$presolve <-
                        as.integer(lpsolver.presolve)
                }
            } else {
                lpsolver.options.criterion <- lpsolver.options.default
            }
            if (hasArg(lpsolver.options.bounds)) {
                lpsolver.options.bounds <-
                    optionsGurobi(lpsolver.options.bounds, debug)
                if (hasArg(lpsolver.presolve)) {
                    lpsolver.options.bounds$presolve <-
                        as.integer(lpsolver.presolve)
                }
            } else {
                lpsolver.options.bounds <- lpsolver.options.default
            }
        }
        ## Turn off non-convex option if using direct regression
        if (direct) lpsolver.options.bounds$nonconvex <- 0
    } else {
        if (lpsolver == "cplexapi") {
            if (hasArg(lpsolver.options)) {
                lpsolver.options.criterion <- optionsCplexAPI(lpsolver.options)
                lpsolver.options.bounds <- optionsCplexAPI(lpsolver.options)
            } else {
                if (hasArg(lpsolver.options.criterion)) {
                    lpsolver.options.criterion <-
                        optionsCplexAPI(lpsolver.options.criterion)
                } else {
                    lpsolver.options.criterion <-
                        optionsCplexAPI(list(setDblParmCPLEX =
                                                 c(parm = 1016,
                                                   value = 1e-06)))
                }
                if (hasArg(lpsolver.options.bounds)) {
                    lpsolver.options.bounds <-
                        optionsCplexAPI(lpsolver.options.bounds)
                } else {
                    lpsolver.options.bounds <-
                        optionsCplexAPI(list(setDblParmCPLEX =
                                                 c(parm = 1016,
                                                   value = 1e-06)))
                }
            }
        }
        if (lpsolver == "lpsolveapi") {
            if (hasArg(lpsolver.options)) {
                lpsolver.options.criterion <-
                    optionsLpSolveAPI(lpsolver.options)
                lpsolver.options.bounds <- optionsLpSolveAPI(lpsolver.options)
            } else {
                if (hasArg(lpsolver.options.criterion)) {
                    lpsolver.options.criterion <-
                        optionsLpSolveAPI(lpsolver.options.criterion)
                } else {
                    lpsolver.options.criterion <-
                        optionsLpSolveAPI(list(epslevel = "tight"))
                }
                if (hasArg(lpsolver.options.bounds)) {
                    lpsolver.options.bounds <-
                        optionsLpSolveAPI(lpsolver.options.bounds)
                } else {
                    lpsolver.options.bounds <-
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
        ## Generate LP environment that is to be updated
        lpEnv <- new.env()
        lpEnv$mbobj <- eval(monoboundAcall)
    }
    ## Setup LP problem
    lpSetup(env = lpEnv, sset = sset, orig.sset = NULL,
            lpsolver = lpsolver, direct = direct)
    ## Setup QCQP problem
    if (direct) qpSetup(env = lpEnv, sset = sset, g0 = gstar0, g1 = gstar1,
                        criterion.tol = criterion.tol,
                        qpsolver = lpsolver)
    ## Prepare LP messages
    ##
    ## Status codes: 0-unknown; 1-optimal; 2-infeasible; 3-infeasible or
    ## unbounded; 4-unbounded; 5-numerical error; 6-suboptimal;
    ## 7-optimal but infeasible after rescaling.
    messageAlt <- gsub("\\s+", " ",
                       "If the LP solver does not return a solution,
                        consider exporting the LP model
                        and passing it to the LP solver
                        for more details.\n")
    messageInf <- gsub("\\s+", " ",
                       "Since a minimum criterion was found, the
                        model should be feasible.
                        For more details, consider exporting the model
                        and passing it to the LP solver.\n")
    messageInfUnb <- gsub("\\s+", " ",
                          "Since a minimum criterion was found, the
                           model is most likely feasible but
                           unbounded. This can happen if
                           the initial grid is too small. Try
                           increasing the parameters 'initgrid.nx'
                           and 'initgrid.nu'. If the model is
                           indeed infeasible,
                           consider exporting the model and passing it
                           to the LP solver for more details.\n")
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
                        LP problems.\n")
    messageSub <- gsub("\\s+", " ",
                       "Tolerance parameters for the LP solver
                        can be passed through the argument
                        'lpsolver.options'.\n")
    messageOptInf <- gsub("\\s+", " ",
                          "A possible reason for this is that covariates
                           are not scaled appropriately
                           (i.e. the range of magnitudes exceeds 1e13).
                           Tolerance parameters for the LP solver
                           can also be passed through the argument
                           'lpsolver.options'.\n")
    while (audit_count <= audit.max) {
        if (noisy) {
            cat("\n    Audit count: ", audit_count, "\n", sep = "")
        }
        if (!direct) {
            lpSetupSolver(env = lpEnv, lpsolver = lpsolver)
            lpSetupCriterion(env = lpEnv, sset = sset)
            minobseq <- criterionMin(lpEnv, sset, lpsolver,
                                     lpsolver.options.criterion, debug)
            ## Try to diagnose cases where the solution is not
            ## available. This could be due to infeasibility or numerical
            ## issues. To deal with infeasibilty, the LP problem is solved
            ## without any shape restrictions. We then check if any of the
            ## lower and upper bounds are violated, which is a likely
            ## cause for infeasible solutions.
            if (minobseq$status %in% c(0, 2, 3, 4, 5)) {
                origMinStatus <- minobseq$status
                ## Stop if issues are numerical, or unbounded, or unknown.
                if (origMinStatus == 0) {
                    errMess <-
                        paste('No solution provided by the LP solver.',
                              messageAlt)
                    stop(errMess, call. = FALSE)
                }
                if (origMinStatus == 4) {
                    errMess <-
                        paste0('No solution since the model is unbounded.',
                               messageUnb)
                    stop(errMess, call. = FALSE)
                }
                if (origMinStatus == 5) {
                    errMess <-
                        paste('No solution due to numerical issues.',
                              messageNum)
                    stop(errMess, call. = FALSE)
                }
                ## Otherwise, continue and test for infeasibility.
                rm(minobseq)
                lpSetupInfeasible(lpEnv, sset)
                if (lpsolver == "gurobi") {
                    if (debug & lpsolver.options.criterion$outputflag == 1) {
                        cat("Infeasibility diagnosis optimization statistics:\n")
                        cat("------------------------------------------------\n")
                    }
                }
                minobseqAlt <- criterionMin(env = lpEnv,
                                            sset = sset,
                                            lpsolver = lpsolver,
                                            lpsolver.options =
                                                lpsolver.options.criterion)
                solVec <- minobseqAlt$x
                ## Test for violation
                negatepos <- which(lpEnv$mbobj$mbs == ">=")
                violateDiff <- lpEnv$mbobj$mbA %*% solVec - lpEnv$mbobj$mbrhs
                violateDiff[negatepos] <- -violateDiff[negatepos]
                violatevec <- as.vector(violateDiff > audit.tol)
                violatepos <- which(violatevec == TRUE)
                violateType <- sapply(violatepos, function(x) {
                    if (x %in% lpEnv$mbobj$lb0seq) {
                        if (m0.lb.default == TRUE) {
                            return(paste0("m0.lb = ", round(m0.lb, 6),
                                          " (min. observed outcome by default)"))
                        } else {
                            return(paste0("m0.lb = ", round(m0.lb, 6)))
                        }
                    }
                    if (x %in% lpEnv$mbobj$lb1seq) {
                        if (m1.lb.default == TRUE) {
                            return(paste0("m1.lb = ", round(m1.lb, 6),
                                          " (min. observed outcome by default)"))
                        } else {
                            return(paste0("m1.lb = ", round(m1.lb, 6)))
                        }
                    }
                    if (x %in% lpEnv$mbobj$lbteseq) {
                        if (mte.lb.default == TRUE) {
                            return(paste0("mte.lb = ", round(mte.lb, 6),
                                          " (min. treatment effect 6by default)"))
                        } else {
                            return(paste0("mte.lb = ", round(mte.lb, 6)))
                        }
                    }
                    if (x %in% lpEnv$mbobj$ub0seq) {
                        if (m0.ub.default == TRUE) {
                            return(paste0("m0.ub = ", round(m0.ub, 6),
                                          " (max. observed outcome by default)"))
                        } else {
                            return(paste0("m0.ub = ", round(m0.ub, 6)))
                        }
                    }
                    if (x %in% lpEnv$mbobj$ub1seq) {
                        if (m1.ub.default == TRUE) {
                            return(paste0("m1.ub = ", round(m1.ub, 6),
                                          " (max. observed outcome by default)"))
                        } else {
                            return(paste0("m1.ub = ", round(m1.ub, 6)))
                        }
                    }
                    if (x %in% lpEnv$mbobj$ubteseq) {
                        if (mte.ub.default == TRUE) {
                            return(paste0("mte.ub = ", round(mte.ub, 6),
                                          " (max. treatment effect by default)"))
                        } else {
                            return(paste0("mte.ub = ", round(mte.ub, 6)))
                        }
                    }
                    if (x %in% lpEnv$mbobj$mono0seq) {
                        if (m0.inc == TRUE) return("m0.inc = TRUE")
                        if (m0.dec == TRUE) return("m0.dec = TRUE")
                    }
                    if (x %in% lpEnv$mbobj$mono1seq) {
                        if (m1.inc == TRUE) return("m1.inc = TRUE")
                        if (m1.dec == TRUE) return("m1.dec = TRUE")
                    }
                    if (x %in% lpEnv$mbobj$monomteseq) {
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
                              paste("No solution since the LP solver proved the
                                 model was infeasible.",
                                 messageInfDiag)),
                         call. = FALSE)
                }
                if (origMinStatus == 3) {
                    stop(gsub("\\s+", " ",
                              paste("No solution since the LP solver proved the
                                 model was infeasible or unbounded.",
                                 messageInfDiag, messageUnbDiag)),
                         call. = FALSE)
                }
            }
            if (noisy) {
                cat("    Minimum criterion: ", fmtResult(minobseq$obj), "\n",
                    sep = "")
            }
            ## Perform specification test
            if (!is.null(orig.sset) & !is.null(orig.criterion)) {
                lpSetupCriterionBoot(lpEnv, sset, orig.sset,
                                     orig.criterion, criterion.tol,
                                     setup = TRUE)
                minobseqTest <- criterionMin(lpEnv, sset, lpsolver,
                                             lpsolver.options.criterion)
                lpSetupCriterionBoot(lpEnv, sset, orig.sset,
                                     orig.criterion, criterion.tol,
                                     setup = FALSE)
            }
            ## Provide warnings if solutions are suboptimal.
            bWarn <- NULL
            if (minobseq$status == 6) {
                bWarn <-
                    paste(bWarn,
                          gsub("\\s+", " ",
                               paste('The LP solver was unable to satisfy
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
            if (!is.null(bWarn)) warning(bWarn, call. = FALSE, immediate. = TRUE)
        }

        ## Obtain bounds
        if (noisy) {
            cat("    Obtaining bounds...\n")
        }
        if (!direct) {
            lpSetupBound(env = lpEnv,
                         g0 = gstar0,
                         g1 = gstar1,
                         sset = sset,
                         criterion.factor = minobseq$obj * (1 + criterion.tol),
                         lpsolver = lpsolver,
                         setup = TRUE)
        }
        lpresult <- bound(env = lpEnv,
                          sset = sset,
                          lpsolver = lpsolver,
                          lpsolver.options = lpsolver.options.bounds,
                          noisy = noisy,
                          smallreturnlist = smallreturnlist,
                          debug = debug)
        if (lpresult$error == TRUE) {
            errMess <- NULL
            errTypes <- NULL
            for (type in c('min', 'max')) {
                tmpName <- paste0(type, 'status')
                if (type == 'min') tmpType <- 'minimization'
                if (type == 'max') tmpType <- 'maximization'
                if (lpresult[[tmpName]] == 0) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The LP solver did not provide a
                                     solution for the', tmpType,
                                     'problem.')))
                }
                if (lpresult[[tmpName]] == 2) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType, 'problem was proven
                                          to be infeasible.')))
                }
                if (lpresult[[tmpName]] == 3) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType,
                                         'problem was proven to be
                                          infeasible or unbounded.')))
                }
                if (lpresult[[tmpName]] == 4) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType,
                                         'problem was proven to be
                                          unbounded.')))
                }
                if (lpresult[[tmpName]] == 5) {
                    errMess <-
                        paste(errMess,
                              gsub('\\s+', ' ',
                                   paste('The', tmpType,
                                         'problem resulted in
                                          numerical issues.')))
                }
                errTypes <- c(errTypes, lpresult[[tmpName]])
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
                        audit.grid = audit.grid))
        }
        ## Provide warnings if solutions are suboptimal.
        bWarn <- NULL
        bWarnTypes <- NULL
        for (type in c('min', 'max')) {
            tmpName <- paste0(type, 'status')
            if (type == 'min') tmpType <- 'minimization'
            if (type == 'max') tmpType <- 'maximization'
            if (lpresult[[tmpName]] == 6) {
                bWarn <-
                    paste(bWarn,
                          gsub("\\s+", " ",
                               paste('The LP solver was unable to satisfy
                               the optimality tolerance for the',
                               tmpType, 'problem, so a suboptimal
                               solution is returned.')))
            }
            if (lpresult[[tmpName]] == 7) {
                bWarn <-
                    paste(bWarn,
                          gsub("\\s+", " ",
                               paste('The solution to the',
                                     tmpType, 'problem is optimal,
                                     but infeasible after rescaling.')))
            }
            bWarnTypes <- c(bWarnTypes,
                            lpresult[[tmpName]])
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
        solVecMin <- c(lpresult$ming0, lpresult$ming1)
        solVecMax <- c(lpresult$maxg0, lpresult$maxg1)
        optstatus <- min(c(lpresult$minstatus,
                           lpresult$maxstatus))
        if (existsolution == FALSE) existsolution <- TRUE
        prevbound <- c(lpresult$min, lpresult$max)

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
                                                     lpresult$ming0,
                                                 solution.m1.min =
                                                     lpresult$ming1,
                                                 solution.m0.max =
                                                     lpresult$maxg0,
                                                 solution.m1.max =
                                                     lpresult$maxg1,
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
                                           the LP constraint grid and
                                           audit grid are identical,
                                           and the audit tolerance has
                                           increased from ",
                                           format(origTol, scientific = TRUE),
                                           " to ",
                                           format(audit.tol, scientific = TRUE),
                                           ". This suggests precision of the
                                           LP solver exceeds that of R,
                                           and may be due to scaling issues
                                           in the LP model.
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
                                         the LP constraint grid and audit grid
                                         being identical.
                                    Audit tolerance was increased from ",
                                    format(origTol, scientific = TRUE), " to ",
                                    format(audit.tol, scientific = TRUE),
                                    ".  This suggests precision of the
                                 LP solver exceeds that of R,
                                 and may be due to scaling issues
                                 in the LP model.")), "\n",
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
                nShapeConstraints <- max(c(lpEnv$mbobj$lb0seq,
                                           lpEnv$mbobj$lb1seq,
                                           lpEnv$mbobj$ub0seq,
                                           lpEnv$mbobj$ub1seq,
                                           lpEnv$mbobj$lbteseq,
                                           lpEnv$mbobj$ubteseq,
                                           lpEnv$mbobj$mono0seq[, 1],
                                           lpEnv$mbobj$mono1seq[, 1],
                                           lpEnv$mbobj$monoteseq[, 1]))
                for (i in types) {
                    ## Expand constraints for m0
                    if (i == 1) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
                        addlb0seq <- seq(length(addIndex)) + tmpAdd
                        tmpAdd <- tmpAdd + length(addIndex)
                        tmpGrid <- auditObj$bounds$bdA$m0.lb
                        lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                               rep(">=", length(addIndex)))
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
                                             rep(m0.lb, length(addIndex)))
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
                        lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                               rep("<=", length(addIndex)))
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
                                             rep(m0.ub, length(addIndex)))
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
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
                            addmono0decseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep("<=", length(addIndex)))
                        }
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
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
                    lpEnv$lpobj$A <-
                        rbind(lpEnv$lpobj$A,
                              cbind(matrix(0, nrow = nrow(addm0),
                                           ncol = 2 * sn),
                                    addm0,
                                    matrix(0, nrow = nrow(addm0),
                                           ncol = addCol)))
                    rm(addCol)
                    ## Update the contraint sequences
                    lpEnv$mbobj$lb0seq <- c(lpEnv$mbobj$lb0seq,
                                            addlb0seq + nShapeConstraints)
                    rm(addlb0seq)
                    lpEnv$mbobj$ub0seq <- c(lpEnv$mbobj$ub0seq,
                                            addub0seq + nShapeConstraints)
                    rm(addub0seq)
                    if (!is.null(addmono0incseq)) {
                        lpEnv$mbobj$mono0seq <-
                            rbind(lpEnv$mbobj$mono0seq,
                                  cbind(addmono0incseq + nShapeConstraints, 1))
                    }
                    if (!is.null(addmono0decseq)) {
                        lpEnv$mbobj$mono0seq <-
                            rbind(lpEnv$mbobj$mono0seq,
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
                        lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                               rep(">=", length(addIndex)))
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
                                             rep(m1.lb, length(addIndex)))
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
                        lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                               rep("<=", length(addIndex)))
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
                                             rep(m1.ub, length(addIndex)))
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
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
                            addmono1decseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep("<=", length(addIndex)))
                        }
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
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
                    lpEnv$lpobj$A <-
                        rbind(lpEnv$lpobj$A,
                              cbind(matrix(0, nrow = nrow(addm1),
                                           ncol = 2 * sn),
                                    matrix(0, nrow = nrow(addm1),
                                           ncol = addCol),
                                    addm1))
                    rm(addCol)
                    ## Update the contraint sequences
                    lpEnv$mbobj$lb1seq <- c(lpEnv$mbobj$lb1seq,
                                            addlb1seq + nShapeConstraints)
                    rm(addlb1seq)
                    lpEnv$mbobj$ub1seq <- c(lpEnv$mbobj$ub1seq,
                                            addub1seq + nShapeConstraints)
                    rm(addub1seq)
                    if (!is.null(addmono1incseq)) {
                        lpEnv$mbobj$mono1seq <-
                            rbind(lpEnv$mbobj$mono1seq,
                                  cbind(addmono1incseq + nShapeConstraints, 1))
                    }
                    if (!is.null(addmono1decseq)) {
                        lpEnv$mbobj$mono1seq <-
                            rbind(lpEnv$mbobj$mono1seq,
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
                        lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                               rep(">=", length(addIndex)))
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
                                             rep(mte.lb, length(addIndex)))
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
                        lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                               rep("<=", length(addIndex)))
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
                                             rep(mte.ub, length(addIndex)))
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
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
                            addmonotedecseq <- seq(length(addIndex)) + tmpAdd
                            tmpAdd <- tmpAdd + length(addIndex)
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep("<=", length(addIndex)))
                        }
                        lpEnv$lpobj$rhs <- c(lpEnv$lpobj$rhs,
                                             rep(0, length(addIndex)))
                        addmte <- rbind(addmte, tmpGrid[addIndex, ])
                        rm(tmpGrid)
                        if (i == 12) auditObj$mono$monoA$mte <- NULL
                    }
                }
                if (!is.null(addmte)) {
                    ## Update the constraint matrix
                    lpEnv$lpobj$A <-
                        rbind(lpEnv$lpobj$A,
                              cbind(matrix(0, nrow = nrow(addmte),
                                           ncol = 2 * sn),
                                    addmte))
                    ## Update the contraint sequences
                    lpEnv$mbobj$lbteseq <- c(lpEnv$mbobj$lbteseq,
                                             addlbteseq + nShapeConstraints)
                    rm(addlbteseq)
                    lpEnv$mbobj$ubteseq <- c(lpEnv$mbobj$ubteseq,
                                             addubteseq + nShapeConstraints)
                    rm(addubteseq)
                    if (!is.null(addmonoteincseq)) {
                        lpEnv$mbobj$monoteseq <-
                            rbind(lpEnv$mbobj$monoteseq,
                                  cbind(addmonoteincseq + nShapeConstraints, 1))
                    }
                    if (!is.null(addmonotedecseq)) {
                        lpEnv$mbobj$monoteseq <-
                            rbind(lpEnv$mbobj$monoteseq,
                                  cbind(addmonotedecseq + nShapeConstraints,
                                        -1))
                    }
                    rm(addmonoteincseq, addmonotedecseq)
                }
                rm(addmte)
                ## Move on to next iteration of the audit
                audit_count <- audit_count + 1
                lpSetupBound(env = lpEnv, setup = FALSE)
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
    ## Return output
    output <- list(max = lpresult$max,
                   min = lpresult$min,
                   lpresult = lpresult,
                   gridobj = list(audit.grid = audit.grid,
                                  violations = violateMat),
                   auditcount = audit_count)
    if (!direct) output$minobseq <- minobseq$obj
    if (!is.null(orig.sset) && !is.null(orig.criterion)) {
        output$spectest = minobseqTest$obj
    }
    return(output)
}

#' Select points from audit grid to add to the constraint grid
#'
#' This function selects which points from the audit grid should be
#' included into the original grid. Both the constraint grid and audit
#' grid are represented as constraints in an LP problem. This function
#' selects which points in the audit grid (i.e. which rows in the
#' audit constraint matrix) should be added to the constraint grid
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
#'     LP problem (i.e. the points to add to the original grid).
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
