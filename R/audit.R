#' Audit procedure
#'
#' This is the wrapper for running the entire audit procedure. This
#' function sets up the LP problem of minimizing the violation of
#' observational equivalence for the set of IV-like estimands, while
#' satisfying boundedness and monotonicity constraints declared by the
#' user. Rather than enforce boundedness and monotonicity hold across
#' the entire support of covariates and unobservables, this procedure
#' enforces the conditions over a subset of points in a grid. This
#' grid corresponds to the set of values the covariates can take, and a
#' subset of values of the unobservable term. The size of this grid is
#' specified by the user in the function arguments. The procedure then
#' goes on to check whether the constraints are satisfied at points
#' off the grid. Any point where either the boundedness or
#' monotonicity constraints are violated are incorporated into the
#' grid, and the process is repeated until the grid incorporates the
#' entire support of the covariates, or until some a maximum number of
#' iterations is reached.
#'
#' @param audit.grid list, contains the A matrix used in the audit
#'     for the original sample, as well as the RHS vector used in the
#'     audit from the original sample.
#' @param save.grid boolean, set to \code{FALSE} by default. Set to
#'     true if the fine grid from the audit should be saved. This
#'     option is used for inference procedure under partial
#'     identification, which uses the fine grid from the original
#'     sample in all bootstrap resamples.
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
#' @param sset a list containing the point estimates and gamma
#' @param gstar0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param gstar1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param orig.sset list, only used for bootstraps. The list caontains
#'     the gamma moments for each element in the S-set, as well as the
#'     IV-like coefficients.
#' @param orig.criterion numeric, only used for bootstraps. The scalar
#'     corresponds to the minimum observational equivalence criterion
#'     from the original sample.
#'
#' @inheritParams ivmteEstimate
#'
#' @return a list. Included in the list is the minimum violation of
#'     observational equivalence of the set of IV-like estimands, as
#'     well as the list of matrices and vectors associated with
#'     solving the LP problem.
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
#' ## obtained from this object by the function 'genSSet'
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
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
audit <- function(data, uname, m0, m1, splinesobj,
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
                  noisy = TRUE, seed = 12345, debug = FALSE) {
    set.seed(seed)
    call  <- match.call()
    lpsolver <- tolower(lpsolver)
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
                        'mte.dec', 'mte.inc')
    if (audit_count == 1) {
        sn <- length(sset)
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
        uvec <- sort(c(0, 1, round(rhalton(initgrid.nu), 8)))
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
                                                 monov = monov))
        lpEnv <- new.env()
        lpEnv$mbobj <- eval(monoboundAcall)
    }

    ## Generate LP environment that is to be updated
    lpSetup(env = lpEnv, sset = sset, orig.sset = NULL,
            lpsolver = lpsolver)
    while (audit_count <= audit.max) {
        if (noisy) {
            cat("\n    Audit count: ", audit_count, "\n", sep = "")
        }
        lpSetupSolver(env = lpEnv, lpsolver = lpsolver)
        lpSetupCriterion(env = lpEnv, sset = sset)
        minobseq <- obsEqMin(lpEnv, sset, lpsolver,
                             lpsolver.options.criterion, debug)
        ## Try to diagnose cases where the solution is
        ## infeasible. Here, the problem is solved without any shape
        ## restrictions. We then check if any of the lower and upper
        ## bounds are violated, which is a likely cause for infeasible
        ## solutions.
        if (!is.numeric(minobseq$obj) || is.na(minobseq$obj) ||
            (lpsolver == "lpsolveapi" && minobseq$status == 0)) {
            rm(minobseq)
            lpSetupInfeasible(lpEnv, sset)
            minobseqAlt <- obsEqMin(env = lpEnv,
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
            stop(gsub("\\s+", " ",
                      paste0("No feasible solution to the criterion minimization
                              problem. This should only happen if the implied
                              parameter space is empty. The likely cause of an
                              empty parameter space is incoherent shape
                              restrictions. For example, if ",
                             paste(unique(violateType), collapse = ", "),
                             " are all set simultaneously. Try changing the
                              shape constraints on the MTR functions.\n")),
                 call. = FALSE)
        }
        if (noisy) {
            cat("    Minimum criterion: ", fmtResult(minobseq$obj), "\n",
                sep = "")
        }
        ## Perform specification test
        if (!is.null(orig.sset) & !is.null(orig.criterion)) {
            lpSetupCriterionBoot(lpEnv, sset, orig.sset,
                                 orig.criterion, criterion.tol, setup = TRUE)
            minobseqTest <- obsEqMin(lpEnv, sset, lpsolver,
                                     lpsolver.options.criterion)
            lpSetupCriterionBoot(lpEnv, sset, orig.sset,
                                 orig.criterion, criterion.tol, setup = FALSE)
        }

        ## Obtain bounds
        if (noisy) {
            cat("    Obtaining bounds...\n")
        }
        lpSetupBound(env = lpEnv,
                     g0 = gstar0,
                     g1 = gstar1,
                     sset = sset,
                     obseq.factor = minobseq$obj * (1 + criterion.tol),
                     lpsolver = lpsolver,
                     setup = TRUE)
        lpresult <- bound(env = lpEnv,
                          sset = sset,
                          lpsolver = lpsolver,
                          lpsolver.options = lpsolver.options.bounds,
                          noisy = noisy,
                          smallreturnlist = smallreturnlist,
                          debug = debug)
        if (is.null(lpresult)) {
            if (noisy) {
                message("    LP solutions are unbounded.")
            }
            return(list(error = "Failure to maximize/minimize.",
                        audit.grid = audit.grid))
        }
        solVecMin <- c(lpresult$ming0, lpresult$ming1)
        solVecMax <- c(lpresult$maxg0, lpresult$maxg1)
        optstatus <- min(c(lpresult$minstatus,
                           lpresult$maxstatus))
        if (optstatus == 0) {
            if (criterion.tol == 0) {
                stop(gsub("\\s+", " ",
                          "Unable to obtain bounds. Try setting criterion.tol
                          to be greater than 0 to allow for model
                          misspecification, or expanding the initial constraint
                          grid size for imposing the shape restrictions
                          (initgrid.nx, initgrid.nu). \n"))
            } else {
                stop(gsub("\\s+", " ",
                          paste0("Bounds extend to +/- infinity.
                          Consider increasing the size of
                          the initial constraint grid (initgrid.nx,
                          initgrid.nu). \n")))
            }
        } else {
            if (existsolution == FALSE) existsolution <- TRUE
            prevbound <- c(lpresult$min, lpresult$max)
        }
        ## Test for violations for minimization problem
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
                                                 audit.tol = audit.tol))
        auditObj<- eval(monoboundAcall)
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
                ## End audit and return violation matrix if max audit is achieved.
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
                for (i in types) {
                    ## Expand constraints for m0
                    if (i == 1) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
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
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
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
                    lpEnv$lpobj$A <-
                        rbind(lpEnv$lpobj$A,
                              cbind(matrix(0, nrow = nrow(addm0),
                                           ncol = 2 * sn),
                                    addm0,
                                    matrix(0, nrow = nrow(addm0),
                                           ncol = length(sset$s1$g1))))
                }
                rm(addm0)
                ## Expand constraints for m1
                for (i in types) {
                    if (i == 2) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
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
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
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
                    lpEnv$lpobj$A <-
                        rbind(lpEnv$lpobj$A,
                              cbind(matrix(0, nrow = nrow(addm1),
                                           ncol = 2 * sn),
                                    matrix(0, nrow = nrow(addm1),
                                           ncol = length(sset$s1$g0)),
                                    addm1))
                }
                rm(addm1)
                ## Expand constraints for mte
                for (i in types) {
                    if (i == 3) {
                        addIndex <- violateMat[violateMat$type == i,
                                               "pos"]
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
                            lpEnv$lpobj$sense <- c(lpEnv$lpobj$sense,
                                                   rep(">=", length(addIndex)))
                        } else {
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
                    lpEnv$lpobj$A <-
                        rbind(lpEnv$lpobj$A,
                              cbind(matrix(0, nrow = nrow(addmte),
                                           ncol = 2 * sn),
                                    addmte))
                }
                rm(addmte)
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
    output <- list(max = lpresult$max,
                   min = lpresult$min,
                   lpresult = lpresult,
                   minobseq = minobseq$obj,
                   gridobj = list(audit.grid = audit.grid,
                                  violations = violateMat),
                   auditcount = audit_count)
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
#'     violation of a shape constriant.
#' @param audit.add integer, the number of points from the audit grid
#'     to add to the intial for each constraint type. For instance, if
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
#'     m0, and whether the constriant is increasing (+1) or decreasing
#'     (-1).
#' @param mono1seq integer matrix, indicates which rows in the audit
#'     constraint matrix correspond to the monotonicity conditions for
#'     m1, and whether the constriant is increasing (+1) or decreasing
#'     (-1).
#' @param monoteseq integer matrix, indicates which rows in the audit
#'     constraint matrix correspond to the monotonicity conditions for
#'     the treatment effect, and whether the constriant is increasing
#'     (+1) or decreasing (-1).
#' @param mbmap integer vector, indexes the X-value associated with
#'     each row in the audit constraint matrix.
#' @return The audit grid is represented using a set of constraint
#'     matrices. Each point in the audit grid corresponds to a set of
#'     rows in the constraint matrices. The function simply returns
#'     the vector of row numbers for the points from the audit grid
#'     whose corresponding constriants should be added to the original
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
