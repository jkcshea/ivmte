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
#' @param audit.grid list, contains the A A matrix used in the audit
#'     for the original sample, as well as the RHS vector used in the
#'     audit from the original sample.
#' @param save.grid boolean, set to \code{FALSE} by default. Set to
#'     true if the fine grid from the audit should be saved. This
#'     option is used for inference procedure under partial
#'     identification, which uses the fine grid from the original
#'     sample in all bootstrap resamples.
#' @param m1.ub.default boolean, default set to TRUE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m0.ub.default boolean, default set to TRUE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m1.lb.default boolean, default set to TRUE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m0.lb.default boolean, default set to TRUE. Indicator for
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
                  audit.max = 25, audit.grid = NULL,
                  save.grid = FALSE,
                  m1.ub, m0.ub, m1.lb, m0.lb,
                  m1.ub.default = FALSE,
                  m0.ub.default = FALSE,
                  m1.lb.default = FALSE,
                  m0.lb.default = FALSE,
                  mte.ub, mte.lb,
                  m0.dec = FALSE, m0.inc = FALSE,
                  m1.dec = FALSE, m1.inc = FALSE,
                  mte.dec = FALSE, mte.inc = FALSE,
                  sset, gstar0, gstar1,
                  orig.sset = NULL, orig.criterion = NULL,
                  obseq.tol = 0, lpsolver,
                  noisy = TRUE, seed = 12345) {
    set.seed(seed)
    call  <- match.call()
    lpsolver <- tolower(lpsolver)
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
    if (!is.null(audit.grid)) support <- audit.grid$support
    if (is.null(audit.grid))  support <- unique(data[, xvars])
    ## check if support is vector or matrix; it can be a vector if
    ## there is only one X term
    if (is.null(dim(support))) {
        support <- data.frame(support)
        colnames(support) <- xvars
    }
    ## deal with case in which there are no covariates.
    if (length(xvars) == 0) {
        noX <- TRUE
    } else {
        noX <- FALSE
        rownames(support) <- seq(1, nrow(support))
        ## Select first iteration of the grid
        full_index <- seq(1, nrow(support))
        initgrid.nx <- min(initgrid.nx, nrow(support))
        audit.nx <- min(audit.nx, nrow(support))
    }
    ## Begin performing the audit
    prevbound <- c(-Inf, Inf)
    existsolution <- FALSE
    audit_count <- 1
    ## Generate a new grid for the audit
    monoboundAlist <- c('sset', 'gstar0', 'gstar1',
                        'm1.ub', 'm0.ub',
                        'm1.lb', 'm0.lb',
                        'mte.ub', 'mte.lb',
                        'm0.dec', 'm0.inc',
                        'm1.dec', 'm1.inc',
                        'mte.dec', 'mte.inc')
    if (audit_count == 1) {
        if (is.null(audit.grid)) {
            sn <- length(sset)
            a_uvec <- sort(c(0, 1, round(rhalton(audit.nu), 8)))
            if (noX) {
                a_grid <- data.frame(a_uvec)
                colnames(a_grid) <- uname
                a_grid_index <- NULL
            } else {
                ## Generate alternate grid from residual indexes
                if (audit.nx == length(full_index)) {
                    a_grid_index <- full_index
                } else {
                    a_grid_index <- sort(
                        sample(full_index,
                               audit.nx,
                               replace = FALSE,
                               prob = replicate(length(full_index),
                               (1 / length(full_index)))))
                }
            }
            ## Generate all monotonicity and boundedness matrices
            ## for the audit
            if (noisy) cat("    Generating audit grid...\n")
            monoboundAcall <- modcall(call,
                                      newcall = genmonoboundA,
                                      keepargs = monoboundAlist,
                                      newargs = list(m0 = m0,
                                                     m1 = m1,
                                                     uname = uname,
                                                     support = support,
                                                     grid_index =
                                                         a_grid_index,
                                                     uvec = a_uvec,
                                                     splinesobj =
                                                         splinesobj,
                                                     monov = monov))
            a_mbobj <- eval(monoboundAcall)
            a_mbobj$support <- support
        } else {
            sn <- length(sset)
            a_mbobj <- audit.grid
            a_grid_index <- unique(audit.grid$gridobj$map)
            if (length(a_grid_index) == 1 && a_grid_index == 0) {
                a_grid_index <- NULL
                noX <- TRUE
            }
        }
        a_mbA <- a_mbobj$mbA[, (2 * sn + 1):ncol(a_mbobj$mbA)]
        negatepos <- which(a_mbobj$mbs == ">=")
        a_mbA[negatepos, ] <- -a_mbA[negatepos, ]
        a_mbrhs <- a_mbobj$mbrhs
        a_mbrhs[negatepos] <- -a_mbrhs[negatepos]
        ## Generate all monotonicity and boundedness matrices for initial grid
        if (noisy) cat("    Generating initial constraint grid...")
        if (noX) {
            grid_index <- NULL
        } else {
            grid_index <- sort(
                sample(a_grid_index,
                       initgrid.nx,
                       replace = FALSE,
                       prob = replicate(length(a_grid_index),
                       (1/length(a_grid_index)))))
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
        mbobj <- eval(monoboundAcall)
    }

    while (audit_count <= audit.max) {
        if (noisy) {
            cat("\n\n    Audit count: ", audit_count, "\n", sep = "")
        }
        ## Minimize violation of observational equivalence
        lpobj <- lpSetup(sset, NULL, mbobj$mbA, mbobj$mbs,
                         mbobj$mbrhs, lpsolver)
        minobseq <- obsEqMin(sset, NULL, NULL,
                             obseq.tol, lpobj, lpsolver)
        ## Try to diagnose cases where the solution is
        ## infeasible. Here, the problem is solved without any shape
        ## restrictions. We then check if any of the lower and upper
        ## bounds are violated, which is a likely cause for infeasible
        ## solutions.
        if (!is.numeric(minobseq$obj) || is.na(minobseq$obj) ||
            (lpsolver == "lpsolveapi" && minobseq$status == 0)) {
            lpobjAlt <- lpSetup(sset = sset,
                                orig.sset = NULL,
                                mbA = mbobj$mbA,
                                mbs = mbobj$mbs,
                                mbrhs = mbobj$mbrhs,
                                lpsolver = lpsolver,
                                shape = FALSE)
            minobseqAlt <- obsEqMin(sset = sset,
                                    orig.sset = NULL,
                                    orig.criterion = NULL,
                                    obseq.tol = obseq.tol,
                                    lpobj = lpobjAlt,
                                    lpsolver = lpsolver)
            solVec <- minobseqAlt$result$x
            ## Test for violations
            mbA <- mbobj$mbA
            negatepos <- which(mbobj$mbs == ">=")
            mbA[negatepos, ] <- -mbA[negatepos, ]
            mbrhs <- mbobj$mbrhs
            mbrhs[negatepos] <- -mbrhs[negatepos]
            violateDiff <- round(mbA %*% solVec - mbrhs,
                                 digits = 10)
            violatevec <- violateDiff > 0
            violatepos <- which(violatevec == TRUE)
            violateType <- sapply(violatepos, function(x) {
                if (x %in% mbobj$lb0seq) {
                    if (m0.lb.default == TRUE) {
                        return(paste0("m0.lb = ", m0.lb,
                                      " (min. observed outcome by default)"))
                    } else {
                        return(paste0("m0.lb = ", m0.lb))
                    }
                }
                if (x %in% mbobj$lb1seq) {
                    if (m1.lb.default == TRUE) {
                        return(paste0("m1.lb = ", m1.lb,
                                      " (min. observed outcome by default)"))
                    } else {
                        return(paste0("m1.lb = ", m1.lb))
                    }
                }
                if (x %in% mbobj$lbteseq) {
                    return(paste0("mte.lb = ", mte.lb))
                }
                if (x %in% mbobj$ub0seq) {
                    if (m0.ub.default == TRUE) {
                        return(paste0("m0.ub = ", m0.ub,
                                      " (max. observed outcome by default)"))
                    } else {
                        return(paste0("m0.ub = ", m0.ub))
                    }
                }
                if (x %in% mbobj$ub1seq) {
                    if (m1.ub.default == TRUE) {
                        return(paste0("m1.ub = ", m1.ub,
                                      " (max. observed outcome by default)"))
                    } else {
                        return(paste0("m1.ub = ", m1.ub))
                    }
                }
                if (x %in% mbobj$ubteseq) {
                    return(paste0("mte.ub = ", mte.ub))
                }
                if (x %in% mbobj$mono0seq) {
                    if (m0.inc == TRUE) return("m0.inc = TRUE")
                    if (m0.dec == TRUE) return("m0.dec = TRUE")
                }
                if (x %in% mbobj$mono1seq) {
                    if (m1.inc == TRUE) return("m1.inc = TRUE")
                    if (m1.dec == TRUE) return("m1.dec = TRUE")
                }
                if (x %in% mbobj$monomteseq) {
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
            lpobjTest <- lpSetup(sset, orig.sset, mbobj$mbA, mbobj$mbs,
                                 mbobj$mbrhs, lpsolver)
            minobseqTest <- obsEqMin(sset, orig.sset, orig.criterion,
                                     obseq.tol, lpobjTest, lpsolver)
        }

        ## Obtain bounds
        if (noisy) {
            cat("    Obtaining bounds...\n")
        }
        lpresult  <- bound(g0 = gstar0,
                           g1 = gstar1,
                           sset = sset,
                           lpobj = lpobj,
                           obseq.factor = minobseq$obj * (1 + obseq.tol),
                           lpsolver = lpsolver)
        if (is.null(lpresult)) {
            message("    LP solutions are unbounded.")
            cat("    Expanding constraint grid, restarting audit.\n")
            return("Failure to maximize/minimize.")
        }
        solVecMin <- c(lpresult$ming0, lpresult$ming1)
        solVecMax <- c(lpresult$maxg0, lpresult$maxg1)
        optstatus <- min(c(lpresult$minstatus,
                           lpresult$maxstatus))
        if (optstatus == 0) {
            if (obseq.tol == 0) {
                stop(gsub("\\s+", " ",
                          "Unable to obtain bounds. Try setting obseq.tol
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
        violateDiffMin <- round(a_mbA %*% solVecMin - a_mbrhs,
                                digits = 10)
        violatevecMin <- violateDiffMin > 0
        ## Test for violations for maximization problem
        violateDiffMax <- round(a_mbA %*% solVecMax - a_mbrhs,
                                digits = 10)
        violatevecMax <- violateDiffMax > 0
        ## Generate violation data set
        violatevec <- violatevecMin + violatevecMax
        violate <- as.logical(sum(violatevec))
        violatevec <- as.logical(violatevec)
        if (violate) {
            if (noisy) cat("    Violations: ", sum(violatevec), "\n")
            ## Store all points that violate the constraints
            diffVec <- mapply(max,
                              violateDiffMin * violatevecMin,
                              violateDiffMax * violatevecMax)
            violateMat <- selectViolations(diffVec = diffVec,
                                           audit.add = audit.add,
                                           lb0seq = a_mbobj$lb0seq,
                                           ub0seq = a_mbobj$ub0seq,
                                           lb1seq = a_mbobj$lb1seq,
                                           ub1seq = a_mbobj$ub1seq,
                                           lbteseq = a_mbobj$lbteseq,
                                           ubteseq = a_mbobj$ubteseq,
                                           mono0seq = a_mbobj$mono0seq,
                                           mono1seq = a_mbobj$mono1seq,
                                           monoteseq = a_mbobj$monoteseq,
                                           mbmap = a_mbobj$mbmap)
            violateIndexes <- violateMat$row
            ## Expand constraint grid
            mbobj$mbA <- rbind(mbobj$mbA, a_mbobj$mbA[violateIndexes, ])
            mbobj$mbrhs <- c(mbobj$mbrhs, a_mbobj$mbrhs[violateIndexes])
            mbobj$mbs <- c(mbobj$mbs, a_mbobj$mbs[violateIndexes])
            audit_count <- audit_count + 1
            if (audit_count <= audit.max) {
                if (noisy) {
                    if (length(violateIndexes) > 1) ps <- 'points'
                    if (length(violateIndexes) == 1) ps <- 'point'
                    cat("    ",
                        gsub("\\s+", " ",
                             paste0("Expanding constraint grid to
                                        include ", length(violateIndexes),
                                    " additional ", ps, "...")), sep = "")
                }
            } else {
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
            }
        } else {
            violateIndexes <- NULL
            if (noisy) {
                cat("    Violations: 0\n")
                cat("    Audit finished.\n\n")
            }
            break
        }
    }
    if (audit_count > audit.max) audit_count <- audit_count - 1
    if (length(violateIndexes) == 0) {
        violations <- NULL
    } else {
        if (noX) {
            violations <- NULL
        } else {
            violations <- support[a_mbobj$mbmap[violateIndexes], ]
        }
        violations <- data.frame(cbind(violations,
                                       matrix(a_mbobj$mbumap[violateIndexes, ],
                                              ncol = 2)[, 2]))
        colnames(violations)[ncol(violations)] <- uname
        violations$.violation.type <- ""
        typeStr <- rep("", nrow(violateMat))
        typeStr[violateMat$type == 1] <- "m0.lb"
        typeStr[violateMat$type == 2] <- "m1.lb"
        typeStr[violateMat$type == 3] <- "mte.lb"
        typeStr[violateMat$type == 4] <- "m0.ub"
        typeStr[violateMat$type == 5] <- "m1.ub"
        typeStr[violateMat$type == 6] <- "mte.ub"
        typeStr[violateMat$type == 7] <- "m0.inc"
        typeStr[violateMat$type == 8] <- "m0.dec"
        typeStr[violateMat$type == 9] <- "m1.inc"
        typeStr[violateMat$type == 10] <- "m1.dec"
        typeStr[violateMat$type == 11] <- "mte.inc"
        typeStr[violateMat$type == 12] <- "mte.dec"
        violations$.violation.type <- typeStr
        rownames(violations) <- seq(nrow(violations))
    }
    output <- list(max = lpresult$max,
                   min = lpresult$min,
                   lpresult = lpresult,
                   minobseq = minobseq$obj,
                   gridobj = list(initial = mbobj$gridobj,
                                  audit = a_mbobj$gridobj,
                                  violations = violations),
                   auditcount = audit_count)
    if (!is.null(orig.sset) && !is.null(orig.criterion)) {
        output$spectest = minobseqTest$obj
    }
    if (save.grid) {
        output$gridobj$a_mbobj <- a_mbobj
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
