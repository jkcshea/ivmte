#' Audit procedure
#'
#' This is the wrapper for running the entire audit procedure. This
#' function sets up the LP problem of minimizing the violation of
#' observational equivalance for the set of IV-like estimands, while
#' satisfying boundedness and monotonicity constraints declared by the
#' user. Rather than enforce boundedness and monotonicity hold across
#' the entire support of covariates and unobservables, this procedure
#' enforces the conditions over a subset of points in a grid. This
#' grid coresponds to the set of values the covariates can take, and a
#' subset of values of the unobservable term. The size of this grid is
#' specified by the user in the function arguments. The procedure then
#' goes on to check whether the constraints are satisfied at points
#' off the grid. Any point where either the boundedness or
#' monotonicity constriants are violated are incoporated into the
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
#' @param audit.Nu number of evenly-spread points of the unobservable
#'     u to use to form the grid in the audit procedure.
#' @param audit.Nx number of `evenly' spread points of the covariates
#'     to use to form the grid in the audit procedure.
#' @param audit.add.x number of points to add to the grid for
#'     covariates in each iteration of the audit procedure.
#' @param audit.add.u number of points to add to the grid for the
#'     unobservables in each iteration of the audit procedure.
#' @param audit.max maximum number of iterations in the audit
#'     procedure.
#' @param audit.tol tolerance for determining when to end the audit
#'     procedure.
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
#' @return a list. Included in the list is the minimum violation of
#'     observational equivalence of the set of IV-like estimands, as
#'     well as the list of matrices and vectors associated with
#'     solving the LP problem.
#'
#' @export
audit.mst <- function(data, uname, m0, m1, splinesobj, vars_mtr, terms_mtr,
                      audit.Nu = 20, audit.Nx = 50,
                      audit.add.x = 5, audit.add.u = 3, audit.max = 5,
                      audit.tol = 1e-08, 
                      m1.ub, m0.ub, m1.lb, m0.lb, mte.ub, mte.lb,
                      m0.dec, m0.inc, m1.dec, m1.inc, mte.dec, mte.inc,
                      sset, gstar0, gstar1) {
   
    call  <- match.call()
    audit <- TRUE

    splines <- list(splinesobj[[1]]$splineslist,
                    splinesobj[[2]]$splineslist)

    ## Update MTR formulas to include all terms that interact with
    ## splines
    for (j in terms_mtr) {
        m0 <- update(m0, as.formula(paste("~ . +",
                                          paste(unique(terms_mtr),
                                                collapse = " + "))))
        m1 <- update(m1, as.formula(paste("~ . +",
                                          paste(unique(terms_mtr),
                                                collapse = " + "))))
    }
    
    ## Obtain name of unobservable variable
    if(hasArg(uname)) {
        if(! suppressWarnings(try(class(uname), silent = TRUE) == "character")){
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
    uvec    <- round(seq(0, 1, length.out = audit.Nu), 8)
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
        ## FIX: Tukey distribution does not seem that great
        quantiles <- apply(support, 1, tukeydist, data = data[, xvars])
        quantiles <- cbind(support, quantiles)
        quantiles <- quantiles[order(quantiles[, "quantiles"]), ]
        rownames(quantiles) <- seq(1, nrow(quantiles))
        support <- data.frame(quantiles[, xvars])
        if (dim(support)[2] == 1) colnames(support) <- xvars

        ## Select first iteration of the grid
        full_index <- seq(1, nrow(support))

        audit.Nx <- min(audit.Nx, nrow(support))
        grid_index <- sample(full_index,
                             audit.Nx,
                             replace = FALSE,
                             prob = replicate(nrow(support), (1/nrow(support))))
        grid_resid <- full_index[!(full_index %in% grid_index)]

    }

    ## Begin performing the audit
    minobseqobj <- Inf    
    audit_count <- 1
    while (audit_count <= audit.max & audit == TRUE) {
        cat("Audit count:", audit_count, "\n")
        if (!noX) {
            if (length(grid_resid) == 0) {
                message("Full support now included as grid.")
            }
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
        lpobj <- lpsetup.mst(sset, mbobj$mbA, mbobj$mbs, mbobj$mbrhs)

        minobseq  <- obseqmin.mst(sset, lpobj)

        message(paste("Minimum observational equivalence deviation:",
                      round(minobseq$obj, 6), "\n"))
        
        if (minobseqobj == 0) {
            if (minobseq$obj == 0) {
                audit <- FALSE
                message(gsub("\\s+", " ",
                             "Audit ending: tolerance level for reduction in
                             deviation of observational equivalence
                             reached.\n"))
                break            
            }
        } else if (minobseqobj == Inf) {
            minobseqobj <- minobseq$obj
        } else {
            if (abs(minobseq$obj - minobseqobj)/minobseqobj < audit.tol) {
                audit <- FALSE
                message(gsub("\\s+", " ",
                             "Audit ending: tolerance level for reduction in
                         deviation of observational equivalence reached.\n"))
                break            
            } else {
                minobseqobj <- minobseq$obj
            }
        }
       
        solutionvec <- c(minobseq$g0, minobseq$g1)
        
        ## Generate a new grid for the audit
        a_uvec <- round(runif(audit.add.u), 8)
        
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
            audit.add.x <- min(audit.add.x, length(grid_resid))
            if (audit.add.x == length(grid_resid)) {
                a_grid_index <- grid_resid
            } else {
                a_grid_index <- sample(grid_resid,
                                        audit.add.x,
                                        replace = FALSE,
                                        prob = replicate(nrow(resid_support),
                                        (1/nrow(resid_support))))
            }
            
            if (length(a_grid_index) == 0) {
                a_grid_index <- grid_index
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
        violatevec <- mapply(">", (a_mbA %*% solutionvec), a_mbrhs)
        violate <- as.logical(sum(violatevec))

        
        if (violate) {
            cat("Expanding audit grid...\n")
            violate_pos <- which(violatevec == TRUE)
            violate_index <- unique(a_mbobj$mbmap[violate_pos])
            grid_index <- c(grid_index, violate_index)
            if (!noX) {
                grid_resid <- grid_resid[!grid_resid %in% violate_index]
            }
            uvec <- sort(unique(c(uvec, c(a_mbobj$mbumap[violate_pos, ]))))
            audit_count <- audit_count + 1
        } else {
            audit <- FALSE
            message(gsub("\\s+", " ",
                         "Audit ending: no violations of monotonicity or
                         boundedness restrictions by points chosen off of the
                         grid defining shape restrictions for the LP
                         problem.\n"))
            break            
        }
    }

    return(list(lpobj    = lpobj,
            minobseq = minobseq))
}
