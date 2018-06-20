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
#' @param data \code{data.frame} used to estimate the treatment effects.
#' @param uname variable name for unobservale used in declaring MTRs.
#' @param m0 one-sided formula for marginal treatment response
#'     function for control group.
#' @param m1 one-sided formula for marginal treatment response
#'     function for treated group.
#' @param audit.Nu number of evenly-spread points of the unobservable u to
#'     use to form the grid in the audit procedure.
#' @param audit.Nx number of `evenly' spread points of the covariates to
#'     use to form the grid in the audit procedure.
#' @param audit.add number of points to add to the grid in each
#'     iteration of the audit procedure.
#' @param audit.max maximum number of iterations in the audit
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
audit.mst <- function(data, uname, m0, m1,
                      audit.Nu = 10, audit.Nx = 10, audit.add = 2, audit.max = 5,
                      m1.ub, m0.ub, m1.lb, m0.lb, mte.ub, mte.lb,
                      m0.dec, m0.inc, m1.dec, m1.inc, mte.dec, mte.inc,
                      sset, gstar0, gstar1) {

    call  <- match.call()
    audit <- TRUE

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
    uvec    <- seq(0, 1, length.out = audit.Nu)
    xvars   <- unique(c(all.vars(m0), all.vars(m1)))
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
        grid <- data.frame(uvec)
        colnames(grid) <- uname
        grid_index <- rownames(grid)
        grid_resid <- c()
        gridobj <- list(grid = grid,
                        map  = replicate(audit.Nu, 1))
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
    audit_count <- 1
    while (audit_count <= audit.max & audit == TRUE) {
        cat("Audit count:", audit_count, "\n")
        if (length(grid_resid) == 0) {
            audit <- FALSE
            message("Full support now included as grid.")
        }
        ## generate the first iteration of the grid
        if (!noX) {
            gridobj <- gengrid.mst(grid_index, support, uvec, uname)
        }

        A0 <- design.mst(formula = m0, data = gridobj$grid)$X
        A0 <- A0[, names(gstar0)]

        A1 <- design.mst(formula = m1, data = gridobj$grid)$X
        A1 <- A1[, names(gstar1)]

        ## generate null objects
        bdA     <- NULL
        monoA   <- NULL

        ## generate matrices for imposing bounds on m0 and m1 and
        ## treatment effects
        if (hasArg(m0.lb) | hasArg(m0.ub) |
            hasArg(m1.lb) | hasArg(m1.lb) |
            hasArg(mte.lb) | hasArg(mte.ub)) {

            boundlist  <- c("m0.lb", "m0.ub",
                            "m1.lb", "m1.ub",
                            "mte.lb", "mte.ub")

            boundAcall <- modcall(call,
                                  newcall = genboundA.mst,
                                  keepargs = boundlist,
                                  newargs = list(A0 = quote(A0),
                                                 A1 = quote(A1),
                                                 sset = quote(sset),
                                                 gridobj = quote(gridobj)))
            bdA <- eval(boundAcall)
        }

        ## Prepare to generate matrices for monotonicity constraints
        if (hasArg(m0.inc)  | hasArg(m0.dec) |
            hasArg(m1.inc)  | hasArg(m1.dec) |
            hasArg(mte.inc) | hasArg(mte.dec)) {

            monolist  <- c("m0.dec", "m0.inc",
                           "m1.dec", "m1.inc",
                           "mte.dec", "mte.inc")
            monoAcall <- modcall(call,
                                 newcall = fullgenmonoA.mst,
                                 keepargs = monolist,
                                 newargs = list(A0 = quote(A0),
                                             A1 = quote(A1),
                                             sset = quote(sset),
                                             gridobj = quote(gridobj),
                                             monov = quote(monov),
                                             gstar0 = quote(gstar0),
                                             gstar1 = quote(gstar1)))
            monoA <- eval(monoAcall)
        }

        ## Now generate the full boundedness/monotonicity matrices and vector
        mbobj <- genfullmbA(bdA, monoA)

        ## Minimize violation of observational equivalence
        lpobj <- lpsetup.mst(sset, mbobj$mbA, mbobj$mbs, mbobj$mbrhs)
        minobseq  <- obseqmin.mst(sset, lpobj)
        message(paste("Minimum observational equivalence deviation:",
            round(minobseq$obj, 6), "\n"))

        ## Now perform the audit (which can only happen if there
        ## remain observations in your empirical support that have not
        ## yet been included in the grid)
        if (length(grid_resid) > 0) {
            negatepos <- which(lpobj$mbs == ">=")
            solutionvec <- c(minobseq$g0, minobseq$g1)

            resid_support <- support[grid_resid, ]
            if(is.null(dim(resid_support))) {
                resid_support <- as.matrix(resid_support)
            }

            ## Generate alternate grid from residual indexes
            audit.add <- min(audit.add, length(grid_resid))
            if (audit.add == length(grid_resid)) {
                newgrid_index <- grid_resid
            } else {
                newgrid_index <- sample(grid_resid,
                                        audit.add,
                                        replace = FALSE,
                                        prob = replicate(nrow(resid_support),
                                        (1/nrow(resid_support))))
            }
            a_gridobj <- gengrid.mst(newgrid_index, support, uvec, uname)
            a_grid    <- a_gridobj$grid
            a_map     <- a_gridobj$map

            a_A0 <- design.mst(formula = m0, data = a_gridobj$grid)$X
            a_A1 <- design.mst(formula = m1, data = a_gridobj$grid)$X

            ## Now generate the constraint matrices for the audit grid

            if (hasArg(m0.lb) | hasArg(m0.ub) |
                hasArg(m1.lb) | hasArg(m1.lb) |
                hasArg(mte.lb) | hasArg(mte.ub)) {

                boundlist  <- c("m0.lb", "m0.ub",
                                "m1.lb", "m1.ub",
                                "mte.lb", "mte.ub")

                boundAcall <- modcall(call,
                                      newcall = genboundA.mst,
                                      keepargs = boundlist,
                                      newargs = c(A0 = quote(a_A0),
                                                  A1 = quote(a_A1),
                                                  sset = quote(sset),
                                                  gridobj = quote(a_gridobj)))
                a_bdA <- eval(boundAcall)
            } else {
                a_bdA <-  NULL
            }

            ## Prepare to generate matrices for monotonicity constraints
            if (hasArg(m0.inc)  | hasArg(m0.dec) |
                hasArg(m1.inc)  | hasArg(m1.dec) |
                hasArg(mte.inc) | hasArg(mte.dec)) {

                monolist  <- c("m0.dec", "m0.inc",
                               "m1.dec", "m1.inc",
                               "mte.dec", "mte.inc")
                monoAcall <- modcall(call,
                                     newcall = fullgenmonoA.mst,
                                     keepargs = monolist,
                                     newargs = list(A0 = quote(a_A0),
                                                    A1 = quote(a_A1),
                                                    sset = quote(sset),
                                                    gridobj = quote(a_gridobj),
                                                    monov = quote(monov),
                                                    gstar0 = quote(gstar0),
                                                    gstar1 = quote(gstar1)))
                a_monoA <- eval(monoAcall)
            } else {
                a_monoA <-  NULL
            }

            a_mbobj <- genfullmbA(a_bdA, a_monoA)
            a_mbA <- a_mbobj$mbA[, (2 * sn + 1) : ncol(a_mbobj$mbA)]
            
            negatepos <- which(a_mbobj$mbs == ">=")

            a_mbA[negatepos, ] <- -a_mbA[negatepos, ]
            a_mbrhs <- a_mbobj$mbrhs
            a_mbrhs[negatepos] <- -a_mbrhs[negatepos]

            ## Test for violations
            violatevec <- mapply(">", (a_mbA %*% solutionvec), a_mbrhs)
            violate <- as.logical(sum(violatevec))

            if (violate) {
                ## if there are violations, include the points that violate
                ## the conditions, and update grid resid.
                cat("Expanding audit grid...\n")
                violate_pos <- which(violatevec == TRUE)
                violate_index <- unique(a_mbobj$mbmap[violate_pos])
                grid_index <- c(grid_index, violate_index)
                grid_resid <- grid_resid[!grid_resid %in% violate_index]
                audit_count <- audit_count + 1
            } else {
                audit <- FALSE
            }
        }
    }
    ## return all the matrices you need to solve the LP
    return(list(lpobj    = lpobj,
                minobseq = minobseq))
}
