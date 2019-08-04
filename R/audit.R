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
#' @param data \code{data.frame} used to estimate the treatment
#'     effects.
#' @param uname name declared by user to represent the unobservable
#'     term in the MTRs.
#' @param m0 one-sided formula for marginal treatment response
#'     function for control group. The unobservable term can be
#'     entered in using splines.
#' @param m1 one-sided formula for marginal treatment response
#'     function for treated group.
#' @param splinesobj list of spline components in the MTRs for treated
#'     and control groups. Spline terms are extracted using
#'     \code{\link{removeSplines}}.
#' @param vars_mtr all variables entering into the MTRs for treated
#'     and control groups.
#' @param terms_mtr0 all terms entering into the MTRs for control
#'     group.
#' @param terms_mtr1 all terms entering into the MTRs for treated
#'     group.
#' @param initgrid.nu number of evenly spread points in the interval
#'     [0, 1] of the unobservable u used to form the grid for imposing
#'     shape restrictions on the MTRs.
#' @param initgrid.nx number of evenly spread points of the covariates
#'     to use to form the grid for imposing shape restrictions on the
#'     MTRs.
#' @param audit.nx number of points on the covariates space to audit
#'     in each iteration of the audit procedure.
#' @param audit.nu number of points in the interval [0, 1],
#'     corresponding to the normalized value of the unobservable term,
#'     to audit in each iteration of the audit procedure.
#' @param audit.add maximum number of points to add to the grids for
#'     imposing each kind of shape constraint. So if there are 5
#'     different kinds of shape constraints, there can be at most
#'     \code{audit.add * 5} additional points added to the grid.
#' @param audit.max maximum number of iterations in the audit
#'     procedure.
#' @param audit.tol tolerance for determining when to end the audit
#'     procedure. Namely, if the percentage change in the upper and
#'     lower bounds both fall below \code{audit.tol} between
#'     iterations of the audit, the audit procedure ends.
#' @param m1.ub numeric value for upper bound on MTR for treated
#'     group.
#' @param m0.ub numeric value for upper bound on MTR for control
#'     group.
#' @param m1.lb numeric value for lower bound on MTR for treated
#'     group.
#' @param m0.lb numeric value for lower bound on MTR for control
#'     group.
#' @param m1.ub.default boolean, default set to TRUE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m0.ub.default boolean, default set to TRUE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m1.lb.default boolean, default set to TRUE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param m0.lb.default boolean, default set to TRUE. Indicator for
#'     whether the value assigned was by the user, or set by default.
#' @param mte.ub numeric value for upper bound on treatment effect
#'     paramter of interest.
#' @param mte.lb numeric value for lower bound on treatment effect
#'     paramter of interest.
#' @param m0.dec logical, equal to TRUE if we want MTR for control
#'     group to be weakly monotone decreasing.
#' @param m0.inc logical, equal to TRUE if we want MTR for control
#'     group to be weakly monotone increasing.
#' @param m1.dec logical, equal to TRUE if we want MTR for treated
#'     group to be weakly monotone decreasing.
#' @param m1.inc logical, equal to TRUE if we want MTR for treated
#'     group to be weakly monotone increasing.
#' @param mte.dec logical, equal to TRUE if we want the MTE to be
#'     weakly monotone decreasing.
#' @param mte.inc logical, equal to TRUE if we want the MTE to be
#'     weakly monotone decreasing.
#' @param gstar0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param gstar1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set.
#' @param obseq.tol tolerance level for how much more the solution is
#'     permitted to violate observational equivalence of the IV-like
#'     estimands. The threshold multiplies the violation of the
#'     observational equivalence, i.e. a threshold of 0 corresponds to
#'     the assumption that the model is correctly specified, and that
#'     any violation of observational equivalence is due to
#'     statistical noise.
#' @param lpsolver name of the linear programming package in R used to
#'     obtain the bounds on the treatment effect.
#' @param noisy boolean, set to TRUE by default. If TRUE, then output
#'     throughout the audit procedure is printed.
#' @param seed integer, the seed that determines the random grid in
#'     the audit procedure.
#' @return a list. Included in the list is the minimum violation of
#'     observational equivalence of the set of IV-like estimands, as
#'     well as the list of matrices and vectors associated with
#'     solving the LP problem.
#'
#' @examples
#'
#' set.seed(10L)
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
#'                          uname = u,
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
                  initgrid.nu = 10, initgrid.nx = 20,
                  audit.nx = 2500, audit.nu = 25, audit.add = 100,
                  audit.max = 25, audit.tol = 1e-08,
                  m1.ub, m0.ub, m1.lb, m0.lb,
                  m1.ub.default = FALSE,
                  m0.ub.default = FALSE,
                  m1.lb.default = FALSE,
                  m0.lb.default = FALSE,
                  mte.ub, mte.lb,
                  m0.dec = FALSE, m0.inc = FALSE,
                  m1.dec = FALSE, m1.inc = FALSE,
                  mte.dec = FALSE, mte.inc = FALSE,
                  sset, gstar0, gstar1, obseq.tol = 0.05, lpsolver,
                  noisy = TRUE, seed = 12345) {
    set.seed(seed)
    call  <- match.call()
    lpsolver <- tolower(lpsolver)    
    ## Clean boolean terms
    terms_mtr0 <- parenthBoolean(terms_mtr0)
    terms_mtr1 <- parenthBoolean(terms_mtr1)
    ## Obtain name of unobservable variable
    if(hasArg(uname)) {
        if (!suppressWarnings(try(class(uname), silent = TRUE) ==
                              "character")) {
            uname <- deparse(substitute(uname))
        }
    } else {
        uname <- "u"
    }
    ## Organize variables
    sn      <- length(sset)
    monov   <- uname ## `monov' is a placeholder name for the monotone
                     ## variable. I use this in case I want to
                     ## generalize the monotonciity restrictions to
                     ## other covariates
    uvec    <- round(seq(0, 1, length.out = initgrid.nu), 8)
    xvars   <- unique(vars_mtr)
    xvars   <- xvars[xvars != uname]
    xvars <- xvars[xvars %in% vars_data]
    otherx  <- xvars[xvars != monov]
    support <- unique(data[, xvars])
    ## check if support is vector or matrix; it can be a vector if
    ## there is only one X term
    if (is.null(dim(support))) {
        support <- data.frame(support)
        colnames(support) <- xvars
    }
    ## deal with case in which there are no covariates.
    if (length(xvars) == 0) {
        noX <- TRUE
        grid_index <- NULL
    } else {
        noX <- FALSE
        rownames(support) <- seq(1, nrow(support))
        ## Select first iteration of the grid
        full_index <- seq(1, nrow(support))
        initgrid.nx <- min(initgrid.nx, nrow(support))
        grid_index <- sample(full_index,
                             initgrid.nx,
                             replace = FALSE,
                             prob = replicate(nrow(support), (1/nrow(support))))
        grid_resid <- full_index[!(full_index %in% grid_index)]
    }
    ## Begin performing the audit
    prevbound <- c(-Inf, Inf)
    existsolution <- FALSE
    audit_count <- 1
    while (audit_count <= audit.max) {
        if (noisy) {
            cat("\n\n    Audit count: ", audit_count, "\n", sep = "")
        }
        ## Generate all monotonicity and boundedness matrices for initial grid
        if (audit_count == 1) {
            cat("    Generating initial grid...\n")
            monoboundAlist <- c('sset', 'gstar0', 'gstar1',
                                'm1.ub', 'm0.ub',
                                'm1.lb', 'm0.lb',
                                'mte.ub', 'mte.lb',
                                'm0.dec', 'm0.inc',
                                'm1.dec', 'm1.inc',
                                'mte.dec', 'mte.inc')
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
        ## Minimize violation of observational equivalence
        lpobj <- lpSetup(sset, mbobj$mbA, mbobj$mbs, mbobj$mbrhs, lpsolver)
        minobseq  <- obsEqMin(sset, lpobj, lpsolver)
        ## Try to diagnose cases where the solution is
        ## infeasible. Here, the problem is solved without any shape
        ## restrictions. We then check if any of the lower and upper
        ## bounds are violated, which is a likely cause for infeasible
        ## solutions.
        if (!is.numeric(minobseq$obj) || is.na(minobseq$obj) ||
            (lpsolver == "lpsolve" && minobseq$status == 0) |
            (lpsolver == "lpsolveapi" && minobseq$status == 0)) {
            lpobjAlt <- lpSetup(sset, mbobj$mbA, mbobj$mbs,
                                    mbobj$mbrhs, lpsolver,
                                    shape = FALSE)
            minobseqAlt <- obsEqMin(sset, lpobjAlt, lpsolver)
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
                        return(paste0("m0.lb (set to min. observed outcome, ",
                                      m0.lb, ", by default)"))
                    } else {
                        return(paste0("m0.lb (set to ",
                                      m0.lb, ")"))
                    }
                }
                if (x %in% mbobj$lb1seq) {
                    if (m1.lb.default == TRUE) {
                        return(paste0("m1.lb (set to min. observed outcome, ",
                                      m1.lb, ", by default)"))
                    } else {
                        return(paste0("m1.lb (set to ",
                                      m1.lb, ")"))
                    }
                }
                if (x %in% mbobj$lbteseq) return("mte.lb")
                if (x %in% mbobj$ub0seq) {
                    if (m0.ub.default == TRUE) {
                        return(paste0("m0.ub (set to max. observed outcome, ",
                                      m0.ub, ", by default)"))
                    } else {
                        return(paste0("m0.ub (set to ",
                                      m0.ub, ")"))
                    }
                }
                if (x %in% mbobj$ub1seq) {
                    if (m1.ub.default == TRUE) {
                        return(paste0("m1.ub (set to max. observed outcome, ",
                                      m1.ub, ", by default)"))
                    } else {
                        return(paste0("m1.ub (set to ",
                                      m1.ub, ")"))
                    }
                }
                if (x %in% mbobj$ubteseq) return("mte.ub")
                if (x %in% mbobj$mono0seq) {
                    if (m0.inc == TRUE) return("m0.inc (set to TRUE)")
                    if (m0.dec == TRUE) return("m0.dec (set to TRUE)")
                }
                if (x %in% mbobj$mono1seq) {
                    if (m1.inc == TRUE) return("m1.inc (set to TRUE)")
                    if (m1.dec == TRUE) return("m1.dec (set to TRUE)")
                }
                if (x %in% mbobj$monomteseq) {
                    if (mte.inc == TRUE) return("mte.inc (set to TRUE)")
                    if (mte.dec == TRUE) return("mte.dec (set to TRUE)")
                }
            })
            stop(gsub("\\s+", " ",
                      paste0("No feasible solution to minimizing violation of
                      observational equivalence. The model may be mispecified.
                      Consider altering the specifications for the MTRs.
                      Infeasible specifications include: ",
                      paste(unique(violateType), collapse = ", "), ".\n")))
        }
        if (noisy) {
            cat("    Minimum criterion: ", fmtResult(minobseq$obj), "\n",
                sep = "")
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
        solVecMin <- c(lpresult$ming0, lpresult$ming1)
        solVecMax <- c(lpresult$maxg0, lpresult$maxg1)
        optstatus <- min(c(lpresult$minstatus,
                           lpresult$maxstatus))
        if (optstatus == 0) {
            if (obseq.tol == 0) {
                stop(gsub("\\s+", " ",
                          "Unable to obtain bounds. Try setting obseq.tol
                          to be greater than 0 to allow for model
                          misspecification, or expanding the initial grid
                          size for imposing the shape restrictions
                          (initgrid.nx, initgrid.nu). \n"))
            } else {
                stop(gsub("\\s+", " ",
                          paste0("Bounds extend to +/- infinity.
                          Consider increasing the size of
                          the initial grid (initgrid.nx, initgrid.nu). \n")))
            }
        } else {
            if (existsolution == FALSE) existsolution <- TRUE
            prevbound <- c(lpresult$min, lpresult$max)
        }
        ## Generate a new grid for the audit
        if (audit_count == 1) {
            a_uvec <- sort(c(round(runif(audit.nu), 8), 0, 1))
            if (noX) {
                a_grid <- data.frame(a_uvec)
                colnames(a_grid) <- uname
                a_grid_index <- NULL
            } else {
                resid_support <- support[grid_resid, ]
                if (is.null(dim(resid_support))) {
                    resid_support <- as.matrix(resid_support)
                }

                ## Generate alternate grid from residual indexes
                audit.nx <- min(audit.nx, length(grid_resid))
                if (audit.nx == length(grid_resid)) {
                    a_grid_index <- grid_resid
                } else {
                    a_grid_index <- sample(grid_resid,
                                           audit.nx,
                                           replace = FALSE,
                                           prob = replicate(nrow(resid_support),
                                           (1 / nrow(resid_support))))
                }
                if (length(a_grid_index) == 0) {
                    a_grid_index <- grid_index
                }
            }
            ## Generate all monotonicity and boundedness matrices for the audit
            cat("    Generating audit grid...\n")
            monoboundAcall <- modcall(call,
                                      newcall = genmonoboundA,
                                      keepargs = monoboundAlist,
                                      newargs = list(m0 = m0,
                                                     m1 = m1,
                                                     uname = uname,
                                                     support = support,
                                                     grid_index = a_grid_index,
                                                     uvec = a_uvec,
                                                     splinesobj = splinesobj,
                                                     monov = monov))
            a_mbobj <- eval(monoboundAcall)
            a_mbA <- a_mbobj$mbA[, (2 * sn + 1) : ncol(a_mbobj$mbA)]
            negatepos <- which(a_mbobj$mbs == ">=")
            a_mbA[negatepos, ] <- -a_mbA[negatepos, ]
            a_mbrhs <- a_mbobj$mbrhs
            a_mbrhs[negatepos] <- -a_mbrhs[negatepos]
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
            cat("    Violations: ", sum(violatevec), "\n")
            ## Store all points that violate the constraints
            diffVec <- violateDiffMin * as.integer(violateDiffMin -
                                                   violateDiffMax > 0) +
                violateDiffMax * as.integer(violateDiffMax - violateDiffMin > 0)
            violateIndexes <- selectViolations(diffVec = diffVec,
                                           audit.add = audit.add,
                                           lb0seq = a_mbobj$lb0seq,
                                           lb1seq = a_mbobj$lb1seq,
                                           ub0seq = a_mbobj$ub0seq,
                                           ub1seq = a_mbobj$ub1seq,
                                           mono0seq = a_mbobj$mono0seq,
                                           mono1seq = a_mbobj$mono1seq,
                                           monoteseq = a_mbobj$monoteseq,
                                           mbmap = a_mbobj$mbmap)

            ## Expand initial grid
            mbobj$mbA <- rbind(mbobj$mbA, a_mbobj$mbA[violateIndexes, ])
            mbobj$mbrhs <- c(mbobj$mbrhs, a_mbobj$mbrhs[violateIndexes])
            mbobj$mbs <- c(mbobj$mbs, a_mbobj$mbs[violateIndexes])
            audit_count <- audit_count + 1
            if (audit_count <= audit.max) {
                if (noisy) {
                    cat("    ",
                        gsub("\\s+", " ",
                             paste0("Expanding initial grids to
                                        include ", length(violateIndexes),
                                    " additional points...")), sep = "")
                }
            } else {
                if (noisy) {
                    cat(gsub("\\s+", " ",
                             paste0("Audit finished: maximum number of
                                 audits (audit.max = ", audit.max,
                                 ") reached.\n")))
                }
                break
            }
        } else {
            if (noisy) {
                cat("    Violations: 0\n")
                cat("    Audit finished.\n\n")
            }
            break
        }
    }
    return(list(max = lpresult$max,
                min = lpresult$min,
                lpresult = lpresult,
                minobseq = minobseq$obj,
                gridobj = mbobj$gridobj,
                auditcount = audit_count))
}


#' Select points from audit grid to add to the initial grid
#'
#' This function selects which points from the audit grid should be
#' included into the original grid. Both the initial grid and audit
#' grid are represented as constraints in an LP problem. This function
#' selects which points in the audit grid (i.e. which rows in the
#' audit constraint matrix) should be added to the initial grid
#' (i.e. should be appended to the initial constraint matrix).
#'
#' @param diffVec numeric vector, with a positive value indicating a
#'     violation of a shape constriant.
#' @param audit.add integer, the number of points from the audit grid
#'     to add to the intial for each constraint type. For instance, if
#'     there are 5 different kinds of constraints imposed, and
#'     \code{audit.add = 5}, then up to 30 points may be added to the
#'     initial grid.
#' @param lb0seq integer vector, indicates which rows in the audit
#'     constraint matrix corresponding to the lower bound for m0.
#' @param lb1seq integer vector, indicates which rows in the audit
#'     constraint matrix corresponding to the lower bound for m1.
#' @param ub0seq integer vector, indicates which rows in the audit
#'     constraint matrix corresponding to the upper bound for m0.
#' @param ub1seq integer vector, indicates which rows in the audit
#'     constraint matrix corresponding to the upper bound for m1.
#' @param mono0seq integer vector, indicates which rows in the audit
#'     constraint matrix corresponding to the monotonicity conditions
#'     for m0.
#' @param mono1seq integer vector, indicates which rows in the audit
#'     constraint matrix corresponding to the monotonicity conditions
#'     for m1.
#' @param monoteseq integer vector, indicates which rows in the audit
#'     constraint matrix corresponding to the monotonicity conditions
#'     for the mte.
#' @param mbmap integer vector, indexes the X-value associated with
#'     each row in the audit constraint matrix.
#' @param mbumap numeric vector, indexes the U-value associated with
#'     each row in the audit constraint matrix.
#' @return The audit grid is represented using a set of constraint
#'     matrices. Each point in the audit grid corresponds to a set of
#'     rows in the constraint matrices. The function simply returns
#'     the vector of row numbers for the points from the audit grid
#'     whose corresponding constriants should be added to the original
#'     LP problem (i.e. the points to add to the original grid).
selectViolations <- function(diffVec, audit.add,
                             lb0seq, lb1seq, ub0seq, ub1seq,
                             mono0seq, mono1seq, monoteseq,
                             mbmap) {
    typeVec <- c(rep(1, times = length(lb0seq)),
                 rep(2, times = length(lb1seq)),
                 rep(3, times = length(ub0seq)),
                 rep(4, times = length(ub1seq)),
                 rep(5, sum(times = mono0seq[, 2] > 0)),
                 rep(6, sum(times = mono0seq[, 2] < 0)),
                 rep(7, sum(times = mono1seq[, 2] > 0)),
                 rep(8, sum(times = mono1seq[, 2] < 0)),
                 rep(9, sum(times = monoteseq[, 2] > 0)),
                 rep(10, sum(times = monoteseq[, 2] < 0)))
    ## Store all points that violate the constraints
    violateMat <- data.frame(cbind(c(lb0seq, lb1seq,
                                     ub0seq, ub1seq,
                                     mono0seq[, 1],
                                     mono1seq[, 1],
                                     monoteseq[, 1]),
                                   typeVec,
                                   mbmap))
    colnames(violateMat) <- c("row", "type", "grid.x")
    violateMat$diff <- diffVec
    violateMat <- violateMat[order(violateMat$type,
                                   violateMat$grid.x,
                                   violateMat$diff), ]
    violateMat$i <- seq(1, nrow(violateMat))
    violateMat[violateMat$diff <= 0, "i"] <- 0
    ## For each point in the X-grid, find the U that violats
    ## the constraints the most
    vmaxu <- aggregate(violateMat$i,
                       by = list(violateMat$type,
                                 violateMat$grid.x),
                       FUN = max)
    vmaxu <- vmaxu$x
    vmaxu <- vmaxu[vmaxu > 0]
    violateMat <- violateMat[vmaxu, ]
    ## Select audit.add number of points to add to the initial grid
    violateMat <- violateMat[order(violateMat$type,
                                   -violateMat$diff), ]
    violateMat$counts <- c(unlist(sapply(table(violateMat$type),
                                         function(x) seq(1, x))))
    violateMat <- violateMat[violateMat$counts <= audit.add, ]
    violateIndexes <- violateMat$row
    return(violateIndexes)
}
