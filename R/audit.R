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
#' @param grid.Nu number of evenly spread points in the interval [0,
#'     1] of the unobservable u used to form the grid for imposing
#'     shape restrictions on the MTRs.
#' @param grid.Nx number of evenly spread points of the covariates to
#'     use to form the grid for imposing shape restrictions on the
#'     MTRs.
#' @param audit.Nx number of points on the covariates space to audit
#'     in each iteration of the audit procedure.
#' @param audit.Nu number of points in the interval [0, 1],
#'     corresponding to the normalized value of the unobservable term,
#'     to audit in each iteration of the audit procedure.
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
#'     there is no violation of observational equivalence.
#' @param lpsolver name of the linear programming package in R used to
#'     obtain the bounds on the treatment effect.
#' @return a list. Included in the list is the minimum violation of
#'     observational equivalence of the set of IV-like estimands, as
#'     well as the list of matrices and vectors associated with
#'     solving the LP problem.
#'
#' @export
audit.mst <- function(data, uname, m0, m1, splinesobj,
                      vars_mtr, terms_mtr0, terms_mtr1,
                      grid.nu = 20, grid.nx = 50,
                      audit.nx = 5, audit.nu = 3, audit.max = 5,
                      audit.tol = 1e-08,
                      m1.ub, m0.ub, m1.lb, m0.lb, mte.ub, mte.lb,
                      m0.dec, m0.inc, m1.dec, m1.inc, mte.dec, mte.inc,
                      sset, gstar0, gstar1, obseq.tol = 1, lpsolver) {

    call  <- match.call()

    splines <- list(splinesobj[[1]]$splineslist,
                    splinesobj[[2]]$splineslist)

    ## Update MTR formulas to include all terms that interact with
    ## splines. This is required for generating the matrices to impose
    ## monotoncity and bounds. The terms that interact wtih the
    ## splines, but do not enter into the MTRs on their own, will be
    ## removed in the function genmonoboundA.
    if (!is.null(m0)) {
        m0 <- update(m0, as.formula(paste("~ . +",
                                          paste(unique(terms_mtr0),
                                                collapse = " + "))))
    } else {
        m0 <- as.formula(paste("~",
                               paste(unique(terms_mtr0),
                                     collapse = " + ")))
    }

    if (!is.null(m1)) {
        m1 <- update(m1, as.formula(paste("~ . +",
                                          paste(unique(terms_mtr1),
                                                collapse = " + "))))
    } else {
        m1 <- as.formula(paste("~",
                               paste(unique(terms_mtr1),
                                     collapse = " + ")))
    }

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
    uvec    <- round(seq(0, 1, length.out = grid.nu), 8)
    xvars   <- unique(vars_mtr)

    xvars   <- xvars[xvars != uname]
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

        ## Obtain quantiles for covariates
        ## FIX: Tukey distribution can be improved.
        quantiles <- apply(support, 1, tukeydist, data = data[, xvars])
        quantiles <- cbind(support, quantiles)
        quantiles <- quantiles[order(quantiles[, "quantiles"]), ]
        rownames(quantiles) <- seq(1, nrow(quantiles))
        support <- data.frame(quantiles[, xvars])
        if (dim(support)[2] == 1) colnames(support) <- xvars

        ## Select first iteration of the grid
        full_index <- seq(1, nrow(support))

        grid.nx <- min(grid.nx, nrow(support))
        grid_index <- sample(full_index,
                             grid.nx,
                             replace = FALSE,
                             prob = replicate(nrow(support), (1/nrow(support))))
        grid_resid <- full_index[!(full_index %in% grid_index)]
    }

    ## Begin performing the audit
    prevbound <- c(-Inf, Inf)
    existsolution <- FALSE
    audit_count <- 1
    while (audit_count <= audit.max) {
        if (obseq.tol > 0 ) {
            cat("Audit count:", audit_count, "\n")
            if (!noX) {
                if (length(grid_resid) == 0) {
                    message("Full support of covariates now included as grid.")
                }
            }
        } else {
            message("\nSkipping audit procedure: obseq.tol set to 0. \n")
        }

        ## Generate all monotonicity and boundedness matrices for initial grid
        monoboundAlist <- c('m0', 'm1',
                            'sset', 'gstar0', 'gstar1',
                            'm1.ub', 'm0.ub',
                            'm1.lb', 'm0.lb',
                            'mte.ub', 'mte.lb',
                            'm0.dec', 'm0.inc',
                            'm1.dec', 'm1.inc',
                            'mte.dec', 'mte.inc')

        monoboundAcall <- modcall(call,
                              newcall = genmonoboundA,
                              keepargs = monoboundAlist,
                              newargs = list(uname = uname,
                                             support = support,
                                             grid_index = grid_index,
                                             uvec = uvec,
                                             splines = splines,
                                             monov = monov))
        mbobj <- eval(monoboundAcall)

        ## Minimize violation of observational equivalence
        lpobj <- lpsetup.mst(sset, mbobj$mbA, mbobj$mbs, mbobj$mbrhs, lpsolver)

        ## print("lpobj")
        ## print(lpobj)

        ## stop("end of testing")

        minobseq  <- obseqmin.mst(sset, lpobj, lpsolver)
        
        if (!is.numeric(minobseq$obj) || is.na(minobseq$obj) ||
            (lpsolver == "lpSolve" && minobseq$status == 0) |
            (lpsolver == "lpSolveAPI" && minobseq$status == 0)) {
            stop(gsub("\\s+", " ",
                      "No feasible solution to minimizing violation of
                      observational equivalence. The model may be mispecified.
                      Consider altering the specifications for the MTRs.\n"))
        }

        if (obseq.tol > 0) {
            message(paste("Minimum observational equivalence deviation:",
                          round(minobseq$obj, 6), "\n"))
        }

        ## Obtain bounds
        message("Obtaining bounds...\n")
        lpresult  <- bound.mst(g0 = gstar0,
                               g1 = gstar1,
                               sset = sset,
                               lpobj = lpobj,
                               obseq.tol = minobseq$obj * obseq.tol,
                               lpsolver = lpsolver)

        solVecMin <- c(lpresult$ming0, lpresult$ming1)
        solVecMax <- c(lpresult$maxg0, lpresult$maxg1)

        optstatus <- min(c(lpresult$minstatus,
                           lpresult$maxstatus))

        if (obseq.tol == 0) {
            if (optstatus == 0) {
                message(gsub("\\s+", " ",
                             "Unable to obtain bounds. Try setting obseq.tol
                             to be greater than 0 to allow for violation of
                             observational equivalence, or expanding the grid
                             size for imposing the shape restrictions
                             (grid.nx, grid.nu). \n"))
            }
            warning(gsub("\\s+", " ",
                         "Setting obseq.tol to 0 allows for no violation of
                         observational equivalence. This assumes the model is
                         correctly specified."))
            break
        }

        ## Generate a new grid for the audit
        a_uvec <- round(runif(audit.nu), 8)

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

        ## Add all audit points to the grid if previous optimization failed
        if (optstatus == 0) {
            existsolution <- FALSE

            message(gsub("\\s+", " ",
                         "Bounds extend to +/- infinity. Expanding grid..."))
            message("")

            grid_index <- c(grid_index, a_grid_index)
            uvec <- c(uvec, a_uvec)
            audit_count <- audit_count + 1

            if (audit_count <= audit.max) {
                next
            } else {
                stop(gsub("\\s+", " ",
                          paste0("Estimation terminated: maximum number of
                          audits (audit.max = ", audit.max, ") reached, but
                          bounds extend to +/- infinity. Either increase the
                          number of audits allowed (audit.max), the size of
                          the initial grid (grid.nx, grid.nu), or the
                          expansion of the grid in the audit procedure
                          (audit.nx, audit.nu). \n")))
            }
        } else {
            if (existsolution == FALSE) {
                existsolution <- TRUE
                prevbound <- c(lpresult$min, lpresult$max)
            } else {
                if (((abs(lpresult$min - prevbound[1]) / prevbound[1]) <
                     audit.tol) &
                    ((abs(lpresult$max - prevbound[2]) / prevbound[2]) <
                     audit.tol)) {
                    message(gsub("\\s+", " ",
                                 "Audit ending: change in bounds falls
                                 below tolerance level.\n"))
                    break
                } else {
                    prevbound <- c(lpresult$min, lpresult$max)
                }
            }
        }

        ## Generate all monotonicity and boundedness matrices for the audit
        monoboundAcall <- modcall(call,
                              newcall = genmonoboundA,
                              keepargs = monoboundAlist,
                              newargs = list(uname = uname,
                                             support = support,
                                             grid_index = a_grid_index,
                                             uvec = a_uvec,
                                             splines = splines,
                                             monov = monov))

        a_mbobj <- eval(monoboundAcall)
        a_mbA <- a_mbobj$mbA[, (2 * sn + 1) : ncol(a_mbobj$mbA)]

        negatepos <- which(a_mbobj$mbs == ">=")
        a_mbA[negatepos, ] <- -a_mbA[negatepos, ]
        a_mbrhs <- a_mbobj$mbrhs
        a_mbrhs[negatepos] <- -a_mbrhs[negatepos]

        ## Test for violations
        violatevecMin <- mapply(">", (a_mbA %*% solVecMin), a_mbrhs)
        violatevecMax <- mapply(">", (a_mbA %*% solVecMax), a_mbrhs)

        violatevec <- violatevecMin + violatevecMax
        violate <- as.logical(sum(violatevec))

        if (violate) {

            violate_pos <- which(violatevec == TRUE)
            violate_index <- unique(a_mbobj$mbmap[violate_pos])
            grid_index <- unique(c(grid_index, violate_index))
            if (!noX) {
                grid_resid <- grid_resid[!grid_resid %in% violate_index]
            }

            uvec <- sort(unique(c(uvec, c(a_mbobj$mbumap[violate_pos, ]))))
            audit_count <- audit_count + 1

            if (audit_count <= audit.max) {
                message("Expanding audit grid...\n")
            } else {
                message(gsub("\\s+", " ",
                             paste0("Audit ending: maximum number of audits
                             (audit.max = ", audit.max, ") reached.\n")))
                break
            }
        } else {
            message(gsub("\\s+", " ",
                         "Audit ending: no violations of monotonicity or
                         boundedness restrictions by points chosen off of the
                         grid defining shape restrictions for the LP
                         problem.\n"))
            break
        }
    }

    return(list(max = lpresult$max,
                min = lpresult$min,
                maxresult = lpresult$maxresult,
                minresult = lpresult$minresult,
                solutionMin = solVecMin,
                solutionMax = solVecMax,
                lpresult = lpresult,
                minobseq = minobseq$obj,
                gridobj = mbobj$gridobj))
}
