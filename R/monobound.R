#' Generating the grid for the audit procedure
#'
#' This function takes in a matrix summarizing the support of the
#' covariates, as well as set of points summarizing the support of the
#' unobservable variable. A Cartesian product of the subset of the
#' support of the covariates and the points in the support of the
#' unobservable generates the grid that is used for the audit
#' procedure.
#'
#' @param index a vector whose elements indicate the rows in the
#'     matrix \code{xsupport} to include in the grid.
#' @param xsupport a matrix containing all the unique combinations of
#'     the covariates included in the MTRs.
#' @param usupport a vector of points in the interval [0, 1],
#'     including 0 and 1. The number of points is decided by the
#'     user. The function generates these points using a Halton
#'     sequence.
#' @param uname name declared by user to represent the unobservable
#'     term.
#' @return a list containing the grid used in the audit; a vector
#'     mapping the elements in the support of the covariates to
#'     \code{index}.
gengrid <- function(index, xsupport, usupport, uname) {
    if (length(usupport) == 0) usupport <- 0
    subsupport <- xsupport[index, ]
    if (is.null(dim(subsupport))) {
        subsupport <- data.frame(subsupport)
        colnames(subsupport) <- colnames(xsupport)
    }
    subsupport$.grid.index <- index
    ## generate a record for which rows correspond to which
    ## index---this will be useful for the audit.
    supportrep <- subsupport[rep(seq(1, nrow(subsupport)),
                                 each = length(usupport)), ]
    uvecrep <- rep(usupport, times = length(index))
    grid <- cbind(supportrep, uvecrep, seq(1, length(uvecrep)))
    rownames(grid) <- grid$.grid.order
    map <- grid$.grid.index
    grid$.grid.index <- NULL
    grid$.u.index <- rep(seq(1, length(usupport)), times = nrow(subsupport))
    colnames(grid) <- c(colnames(xsupport), uname, ".grid.order", ".u.order")
    return(list(grid = grid,
                map = map))
}

#' Generating the constraint matrix
#'
#' This function generates the component of the constraint matrix in
#' the LP/QCQP problem pertaining to the lower and upper bounds on the MTRs
#' and MTEs. These bounds are declared by the user.
#' @param A0 the matrix of values from evaluating the MTR for control
#'     observations over the grid generated to perform the audit. This
#'     matrix will be incorporated into the final constraint matrix
#'     for the bounds.
#' @param A1 the matrix of values from evaluating the MTR for control
#'     observations over the grid generated to perform the audit. This
#'     matrix will be incorporated into the final constraint matrix
#'     for the bounds.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set.
#' @param gridobj a list containing the grid over which the
#'     monotonicity and boundedness conditions are imposed on.
#' @param uname name declared by user to represent the unobservable
#'     term.
#' @param m0.lb scalar, lower bound on MTR for control group.
#' @param m0.ub scalar, upper bound on MTR for control group.
#' @param m1.lb scalar, lower bound on MTR for treated group.
#' @param m1.ub scalar, upper bound on MTR for treated group.
#' @param mte.lb scalar, lower bound on MTE.
#' @param mte.ub scalar, upper bound on MTE.
#' @param solution.m0.min vector, the coefficients for the MTR for
#'     \code{D = 0} corresponding to the lower bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m1.min vector, the coefficients for the MTR for
#'     \code{D = 1} corresponding to the lower bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m0.max vector, the coefficients for the MTR for
#'     \code{D = 0} corresponding to the upper bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m1.max vector, the coefficients for the MTR for
#'     \code{D = 1} corresponding to the upper bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param audit.tol feasibility tolerance when performing the
#'     audit. By default to set to be equal \code{1e-06}. This
#'     parameter should only be changed if the feasibility tolerance
#'     of the solver is changed, or if numerical issues result in
#'     discrepancies between the solver's feasibility check and the
#'     audit.
#' @param qp boolean, set to \code{TRUE} if the direct MTR
#'     regression via QP is used.
#' @return a constraint matrix for the LP/QCQP problem, the associated
#'     vector of inequalities, and the RHS vector in the inequality
#'     constraint. The objects pertain only to the boundedness
#'     constraints declared by the user.
genboundA <- function(A0, A1, sset, gridobj, uname, m0.lb, m0.ub,
                      m1.lb, m1.ub, mte.lb, mte.ub,
                      solution.m0.min = NULL, solution.m1.min = NULL,
                      solution.m0.max = NULL, solution.m1.max = NULL,
                      audit.tol, qp = FALSE) {
    if (!is.null(solution.m0.min) && !is.null(solution.m1.min) &&
        !is.null(solution.m0.max) && !is.null(solution.m1.max)) {
        audit <- TRUE
    } else {
        audit <- FALSE
    }
    grid <- gridobj$grid
    gridmap <- gridobj$map
    namesA0 <- colnames(A0)
    namesA1 <- colnames(A1)
    if (!qp) {
        sn <- length(sset)
        namesA  <- c(seq(1, 2 * sn),
                     namesA0,
                     namesA1)
    } else {
        sn <- 0
        namesA <- c(namesA0, namesA1)
    }
    ## Generate place holders for the matrices representing monotonicity
    lbdA0  <- NULL
    lbdA1  <- NULL
    lbdAte <- NULL
    ubdA0  <- NULL
    ubdA1  <- NULL
    ubdAte <- NULL
    m0ub  <- NULL
    m0ubs <- NULL
    m0lb  <- NULL
    m0lbs <- NULL
    m1ub  <- NULL
    m1ubs <- NULL
    m1lb  <- NULL
    m1lbs <- NULL
    telb  <- NULL
    telbs <- NULL
    teub  <- NULL
    teubs <- NULL
    lbdA0seq  <- NULL
    lbdA1seq  <- NULL
    lbdAteseq <- NULL
    ubdA0seq  <- NULL
    ubdA1seq  <- NULL
    ubdAteseq <- NULL
    map  <- NULL
    umap <- NULL
    bdA <- NULL
    if (audit) {
        bdA <- list()
        diff <- NULL
    }
    ## Construct lower bound matrices
    if (hasArg(m0.lb)) {
        duplicatePos <- duplicated(A0)
        uniqueA0 <- matrix(A0[!duplicatePos, ], nrow = sum(!duplicatePos))
        nonDuplicateIndex <- seq(nrow(A0))[!duplicatePos]
        if (!audit) {
            bdA <- rbind(bdA,
                         cbind(matrix(0, nrow = nrow(uniqueA0), ncol = 2 * sn),
                               uniqueA0,
                               matrix(0, nrow = nrow(uniqueA0),
                                      ncol = ncol(A1))))
            m0lb  <- replicate(nrow(uniqueA0), m0.lb)
            m0lbs <- replicate(nrow(uniqueA0), ">=")
            map <- c(map, gridmap[!duplicatePos])
            umap <- c(umap, grid[!duplicatePos, uname])
            lbdA0seq <- seq(1, nrow(uniqueA0))
        } else {
            violateDiff <- cbind(-(uniqueA0 %*% solution.m0.min - m0.lb),
                                 -(uniqueA0 %*% solution.m0.max - m0.lb))
            violateDiff <- apply(violateDiff, 1, max)
            violatePos <- violateDiff > audit.tol
            diff <- c(diff, violateDiff[violatePos])
            bdA$m0.lb <- matrix(uniqueA0[violatePos, ],
                                   nrow = sum(violatePos))
            map <- c(map, gridmap[nonDuplicateIndex[violatePos]])
            umap <- c(umap, grid[nonDuplicateIndex[violatePos], uname])
            if (sum(violatePos) > 0) {
                lbdA0seq <- seq(sum(violatePos))
            } else {
                lbdA0seq <- NULL
            }
        }
        rm(uniqueA0)
    }
    if (hasArg(m1.lb)) {
        duplicatePos <- duplicated(A1)
        uniqueA1 <- matrix(A1[!duplicatePos, ], nrow = sum(!duplicatePos))
        nonDuplicateIndex <- seq(nrow(A1))[!duplicatePos]
        if (!audit) {
            bdA <- rbind(bdA,
                         cbind(matrix(0, nrow = nrow(uniqueA1), ncol = 2 * sn),
                               matrix(0, nrow = nrow(uniqueA1),
                                      ncol = ncol(A0)),
                               uniqueA1))
            m1lb  <- replicate(nrow(uniqueA1), m1.lb)
            m1lbs <- replicate(nrow(uniqueA1), ">=")
            map <- c(map, gridmap[!duplicatePos])
            umap <- c(umap, grid[!duplicatePos, uname])
            lbdA1seq <- seq(1, nrow(uniqueA1))
        } else {
            violateDiff <- cbind(-(uniqueA1 %*% solution.m1.min - m1.lb),
                                 -(uniqueA1 %*% solution.m1.max - m1.lb))
            violateDiff <- apply(violateDiff, 1, max)
            violatePos <- violateDiff > audit.tol
            diff <- c(diff, violateDiff[violatePos])
            bdA$m1.lb <- matrix(uniqueA1[violatePos, ],
                                   nrow = sum(violatePos))
            map <- c(map, gridmap[nonDuplicateIndex[violatePos]])
            umap <- c(umap, grid[nonDuplicateIndex[violatePos], uname])
            if (sum(violatePos) > 0) {
                lbdA1seq <- seq(sum(violatePos))
            } else {
                lbdA1seq <- NULL
            }
        }
        rm(uniqueA1)
    }
    if (hasArg(mte.lb)) {
        if (!audit) {
            bdA <- rbind(bdA,
                         cbind(matrix(0, nrow = nrow(A0), ncol = 2 * sn),
                               -A0, A1))
            telb  <- replicate(nrow(A1), mte.lb)
            telbs <- replicate(nrow(A1), ">=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            lbdAteseq <- seq(1, nrow(A1))
        } else {
            violateDiff <-
                cbind(-(-A0 %*% solution.m0.min +
                        A1 %*% solution.m1.min - mte.lb),
                      -(-A0 %*% solution.m0.max +
                        A1 %*% solution.m1.max - mte.lb))
            violateDiff <- apply(violateDiff, 1, max)
            violatePos <- violateDiff > audit.tol
            diff <- c(diff, violateDiff[violatePos])
            bdA$mte.lb <- matrix(cbind(A0[violatePos, ],
                                          A1[violatePos, ]),
                                    nrow = sum(violatePos))
            map <- c(map, gridmap[violatePos])
            umap <- c(umap, grid[violatePos, uname])
            if (sum(violatePos) > 0) {
                lbdAteseq <- seq(sum(violatePos))
            } else {
                lbdAteseq <- NULL
            }
        }
    }
    ## Construct upper bound matrices
    if (hasArg(m0.ub)) {
        duplicatePos <- duplicated(A0)
        uniqueA0 <- matrix(A0[!duplicatePos, ], nrow = sum(!duplicatePos))
        nonDuplicateIndex <- seq(nrow(A0))[!duplicatePos]
        if (!audit) {
            bdA <- rbind(bdA,
                         cbind(matrix(0, nrow = nrow(uniqueA0), ncol = 2 * sn),
                               uniqueA0,
                               matrix(0, nrow = nrow(uniqueA0),
                                      ncol = ncol(A1))))
            m0ub  <- replicate(nrow(uniqueA0), m0.ub)
            m0ubs <- replicate(nrow(uniqueA0), "<=")
            map <- c(map, gridmap[!duplicatePos])
            umap <- c(umap, grid[!duplicatePos, uname])
            ubdA0seq <- seq(1, nrow(uniqueA0))
        } else {
            violateDiff <- cbind(uniqueA0 %*% solution.m0.min - m0.ub,
                                 uniqueA0 %*% solution.m0.max - m0.ub)
            violateDiff <- apply(violateDiff, 1, max)
            violatePos <- violateDiff > audit.tol
            if (sum(violatePos) > 0) {
                diff <- c(diff, violateDiff[violatePos])
                bdA$m0.ub <- matrix(uniqueA0[violatePos, ],
                                       nrow = sum(violatePos))
                map <- c(map, gridmap[nonDuplicateIndex[violatePos]])
                umap <- c(umap, grid[nonDuplicateIndex[violatePos], uname])
                ubdA0seq <- seq(sum(violatePos))
            } else {
                ubdA0seq <- NULL
            }
        }
        rm(uniqueA0)
    }
    if (hasArg(m1.ub)) {
        duplicatePos <- duplicated(A1)
        uniqueA1 <- matrix(A1[!duplicatePos, ], nrow = sum(!duplicatePos))
        nonDuplicateIndex <- seq(nrow(A1))[!duplicatePos]
        if (!audit) {
            bdA <- rbind(bdA,
                         cbind(matrix(0, nrow = nrow(uniqueA1), ncol = 2 * sn),
                               matrix(0, nrow = nrow(uniqueA1),
                                      ncol = ncol(A0)),
                               uniqueA1))
            m1ub  <- replicate(nrow(uniqueA1), m1.ub)
            m1ubs <- replicate(nrow(uniqueA1), "<=")
            map <- c(map, gridmap[!duplicatePos])
            umap <- c(umap, grid[!duplicatePos, uname])
            ubdA1seq <- seq(1, nrow(uniqueA1))
        } else {
            violateDiff <- cbind(uniqueA1 %*% solution.m1.min - m1.ub,
                                 uniqueA1 %*% solution.m1.max - m1.ub)
            violateDiff <- apply(violateDiff, 1, max)
            violatePos <- violateDiff > audit.tol
            if (sum(violatePos) > 0) {
                diff <- c(diff, violateDiff[violatePos])
                bdA$m1.ub <- matrix(uniqueA1[violatePos, ],
                                       nrow = sum(violatePos))
                map <- c(map, gridmap[nonDuplicateIndex[violatePos]])
                umap <- c(umap, grid[nonDuplicateIndex[violatePos], uname])
                ubdA1seq <- seq(sum(violatePos))
            } else {
                ubdA1seq <- NULL
            }
        }
        rm(uniqueA1)
    }
    if(hasArg(mte.ub)) {
        if (!audit) {
            bdA <- rbind(bdA,
                         cbind(matrix(0, nrow = nrow(A0), ncol = 2 * sn),
                               -A0, A1))
            teub  <- replicate(nrow(A1), mte.ub)
            teubs <- replicate(nrow(A1), "<=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            ubdAteseq <- seq(1, nrow(A1))
        } else {
            violateDiff <-
                cbind(-A0 %*% solution.m0.min + A1 %*% solution.m1.min - mte.ub,
                      -A0 %*% solution.m0.max + A1 %*% solution.m1.max - mte.ub)
            violateDiff <- apply(violateDiff, 1, max)
            violatePos <- violateDiff > audit.tol
            if (sum(violatePos) > 0) {
                diff <- c(diff, violateDiff[violatePos])
                bdA$mte.ub <- matrix(cbind(A0[violatePos, ],
                                              A1[violatePos, ]),
                                        nrow = sum(violatePos))
                map <- c(map, gridmap[violatePos])
                umap <- c(umap, grid[violatePos, uname])
                ubdAteseq <- seq(sum(violatePos))
            } else {
                ubdAteseq <- NULL
            }
        }
    }
    if (!audit) {
        ## Update indexes for types of boundedness constraints
        countseq <- 0
        if (!is.null(lbdA0seq)) {
            countseq <- countseq + length(lbdA0seq)
        }
        if (!is.null(lbdA1seq)) {
            lbdA1seq <- lbdA1seq + countseq
            countseq <- countseq + length(lbdA1seq)
        }
        if (!is.null(lbdAteseq)) {
            lbdAteseq <- lbdAteseq + countseq
            countseq <- countseq + length(lbdAteseq)
        }
        if (!is.null(ubdA0seq)) {
            ubdA0seq <- ubdA0seq + countseq
            countseq <- countseq + length(ubdA0seq)
        }
        if (!is.null(ubdA1seq)) {
            ubdA1seq <- ubdA1seq + countseq
            countseq <- countseq + length(ubdA1seq)
        }
        if (!is.null(ubdAteseq)) {
            ubdAteseq <- ubdAteseq + countseq
            countseq <- countseq + length(ubdAteseq)
        }
        ## Combine remaining vectors and return
        bds   <- c(m0lbs, m1lbs, telbs,
                   m0ubs, m1ubs, teubs)
        bdrhs <- c(m0lb,  m1lb,  telb,
                   m0ub,  m1ub,  teub)
        map <- matrix(map, ncol = 1)
        colnames(map) <- "grid.X.index"
        return(list(A = bdA,
                    sense = bds,
                    rhs = bdrhs,
                    map = map,
                    umap = umap,
                    lb0seq  = lbdA0seq,
                    lb1seq  = lbdA1seq,
                    lbteseq = lbdAteseq,
                    ub0seq  = ubdA0seq,
                    ub1seq  = ubdA1seq,
                    ubteseq = ubdAteseq))
    } else {
        ## Construct violation matrix
        if (sum(unlist(lapply(bdA, function(x) nrow(x)))) > 0) {
            violateMat <- data.frame(pos = c(lbdA0seq,
                                             lbdA1seq,
                                             lbdAteseq,
                                             ubdA0seq,
                                             ubdA1seq,
                                             ubdAteseq),
                                     type = c(rep(1, length(lbdA0seq)),
                                              rep(2, length(lbdA1seq)),
                                              rep(3, length(lbdAteseq)),
                                              rep(4, length(ubdA0seq)),
                                              rep(5, length(ubdA1seq)),
                                              rep(6, length(ubdAteseq))),
                                     grid.x = map,
                                     grid.u = umap,
                                     diff = diff)
            violateMat$group.name <- paste0(violateMat$type, ".",
                                            violateMat$grid.x)
            violateMat$type.string <- factor(violateMat$type,
                                            levels = seq(6),
                                            labels = c('m0.lb',
                                                       'm1.lb',
                                                       'mte.lb',
                                                       'm0.ub',
                                                       'm1.ub',
                                                       'mte.ub'))
            return(list(bdA = bdA,
                        violateMat = violateMat))
        } else {
            return(NULL)
        }
    }
}

#' Generate components of the monotonicity constraints
#'
#' This function generates the matrix and vectors associated with the
#' monotonicity constraints declared by the user. It takes in a grid
#' of the covariates on which the shape constraints are defined, and then
#' calculates the values of the MTR and MTE over the grid. The
#' matrices characterizing the monotonicity conditions can then be
#' obtained by taking first differences over the grid of the
#' unobservable term, within each set of values in the grid of
#' covariate values.
#' @param A0 the matrix of values from evaluating the MTR for control
#'     observations over the grid generated to perform the audit. This
#'     matrix will be incorporated into the final constraint matrix
#'     for the monotonicity conditions.
#' @param A1 the matrix of values from evaluating the MTR for control
#'     observations over the grid generated to perform the audit. This
#'     matrix will be incorporated into the final constraint matrix
#'     for the monotonicity conditions.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set.
#' @param uname Name of unobserved variable.
#' @param gridobj a list containing the grid over which the
#'     monotonicity and boundedness conditions are imposed on.
#' @param gstar0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param gstar1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param m0.dec boolean, indicating whether the MTR for the control
#'     group is monotone decreasing.
#' @param m0.inc boolean, indicating whether the MTR for the control
#'     group is monotone increasing.
#' @param m1.dec boolean, indicating whether the MTR for the treated
#'     group is monotone decreasing.
#' @param m1.inc boolean, indicating whether the MTR for the treated
#'     group is monotone increasing.
#' @param mte.dec boolean, indicating whether the MTE is monotone
#'     decreasing.
#' @param mte.inc boolean, indicating whether the MTE is monotone
#'     increasing.
#' @param solution.m0.min vector, the coefficients for the MTR for
#'     \code{D = 0} corresponding to the lower bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m1.min vector, the coefficients for the MTR for
#'     \code{D = 1} corresponding to the lower bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m0.max vector, the coefficients for the MTR for
#'     \code{D = 0} corresponding to the upper bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m1.max vector, the coefficients for the MTR for
#'     \code{D = 1} corresponding to the upper bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param audit.tol feasibility tolerance when performing the
#'     audit. By default to set to be equal \code{1e-06}. This
#'     parameter should only be changed if the feasibility tolerance
#'     of the solver is changed, or if numerical issues result in
#'     discrepancies between the solver's feasibility check and the
#'     audit.
#' @param qp boolean, set to \code{TRUE} if the direct MTR
#'     regression via QP is used.
#' @return constraint matrix for the LP/QCQP problem. The matrix pertains
#'     only to the monotonicity conditions on the MTR and MTE declared
#'     by the user.
genmonoA <- function(A0, A1, sset, uname, gridobj, gstar0, gstar1,
                     m0.dec, m0.inc, m1.dec, m1.inc, mte.dec,
                     mte.inc,
                     solution.m0.min = NULL, solution.m1.min = NULL,
                     solution.m0.max = NULL, solution.m1.max = NULL,
                     audit.tol, qp) {
    if (!is.null(solution.m0.min) && !is.null(solution.m1.min) &&
        !is.null(solution.m0.max) && !is.null(solution.m1.max)) {
        audit <- TRUE
    } else {
        audit <- FALSE
    }
    un <- length(unique(gridobj$grid[, uname]))
    grid <- gridobj$grid
    gridmap <- gridobj$map
    ## Construct index for calculating first differences
    uMaxIndex <- unlist(sapply(X = unique(gridmap), FUN = function(x) {
        pos <- sort(which(gridmap == x))
        pos[-1]
    }))
    uMinIndex <- unlist(sapply(X = unique(gridmap), FUN = function(x) {
        pos <- sort(which(gridmap == x))
        pos[-length(pos)]
    }))
    ## Generate list of all relevant matrices. All of these matrices
    ## will be updated along the way.
    monoList <- list(
        monoA0  = NULL,
        monoA1  = NULL,
        monoAte = NULL,
        mono0z  = NULL,
        mono1z  = NULL,
        monotez = NULL,
        mono0s  = NULL,
        mono1s  = NULL,
        monotes = NULL,
        monoA0seq = NULL,
        monoA1seq = NULL,
        monoAteseq = NULL,
        countseq = 0,
        monomap = NULL,
        umap = NULL
    )
    map <- NULL
    umap <- NULL
    ## This matrix should include all the additions 0s on the left
    ## columns
    namesA0 <- colnames(A0)
    namesA1 <- colnames(A1)
    if (!qp) {
        sn <- length(sset)
        namesA  <- c(seq(1, 2 * sn),
                     namesA0,
                     namesA1)
    } else {
        sn <- 0
        namesA  <- c(namesA0,
                     namesA1)
        drY <- sset$s1$ys
        drX <- cbind(sset$s1$g0, sset$s1$g1)
        drQ <- sset$s1$Q
    }
    ## The functions below generate the constraint matrix, the sense
    ## vector, and the RHS vector associated with the monotonicity
    ## constraints for m0, m1, and the mte. In addition, mappings to
    ## the grid index and U index are constructed; also, a list of the
    ## direction of monotonicity constraints is generated (since you
    ## allow m0, m1, and mte to face both increasing and decreasing
    ## monotonicity constraints simultaneously---forcing them to be
    ## constants).
    genmonoA0 <- function(monoObjects, type, audit = FALSE) {
        monoA0 <- A0[uMaxIndex, ] - A0[uMinIndex, ]
        if (is.null(dim(monoA0))) monoA0 <- matrix(monoA0, nrow = 1)
        duplicatePos <- duplicated(monoA0)
        monoA0 <- matrix(monoA0[!duplicatePos, ], nrow = sum(!duplicatePos))
        if (!audit) {
            monoA0 <- cbind(matrix(0, nrow = nrow(monoA0), ncol = 2 * sn),
                            monoA0,
                            matrix(0, nrow = nrow(monoA0), ncol = ncol(A1)))
            colnames(monoA0) <- namesA
        }
        monoObjects$monoA0 <- rbind(monoObjects$monoA0, monoA0)
        monoObjects$mono0z <- c(monoObjects$mono0z, replicate(nrow(monoA0), 0))
        if (type == 1) {
            monoObjects$mono0s <- c(monoObjects$mono0s,
                                    replicate(nrow(monoA0), ">="))
        }
        if (type == -1) {
            monoObjects$mono0s <- c(monoObjects$mono0s,
                                    replicate(nrow(monoA0), "<="))
        }
        monoObjects$monoA0seq <- rbind(monoObjects$monoA0seq,
                                       cbind(seq(1, nrow(monoA0)) +
                                             monoObjects$countseq,
                                             type))
        monoObjects$countseq <- monoObjects$countseq + nrow(monoA0)
        monoObjects$monomap <- c(monoObjects$monomap,
                                 gridmap[uMinIndex[!duplicatePos]])
        monoObjects$umap <- rbind(monoObjects$umap,
                                  cbind(grid[uMinIndex[!duplicatePos], uname],
                                        grid[uMaxIndex[!duplicatePos], uname]))
        return(monoObjects)
    }
    genmonoA1 <- function(monoObjects, type, audit = FALSE) {
        monoA1 <- A1[uMaxIndex, ] - A1[uMinIndex, ]
        if (is.null(dim(monoA1))) monoA1 <- matrix(monoA1, nrow = 1)
        duplicatePos <- duplicated(monoA1)
        monoA1 <- matrix(monoA1[!duplicatePos, ], nrow = sum(!duplicatePos))
        if (!audit) {
            monoA1 <- cbind(matrix(0, nrow = nrow(monoA1), ncol = 2 * sn),
                            matrix(0, nrow = nrow(monoA1), ncol = ncol(A0)),
                            monoA1)
            colnames(monoA1) <- namesA
        }
        monoObjects$monoA1 <- rbind(monoObjects$monoA1, monoA1)
        monoObjects$mono1z <- c(monoObjects$mono1z, replicate(nrow(monoA1), 0))
        if (type == 1) {
            monoObjects$mono1s <- c(monoObjects$mono1s,
                                    replicate(nrow(monoA1), ">="))
        }
        if (type == -1) {
            monoObjects$mono1s <- c(monoObjects$mono1s,
                                    replicate(nrow(monoA1), "<="))
        }
        monoObjects$monoA1seq <- rbind(monoObjects$monoA1seq,
                                       cbind(seq(1, nrow(monoA1)) +
                                             monoObjects$countseq,
                                             type))
        monoObjects$countseq <- monoObjects$countseq + nrow(monoA1)
        monoObjects$monomap <- c(monoObjects$monomap,
                                 gridmap[uMinIndex[!duplicatePos]])
        monoObjects$umap <- rbind(monoObjects$umap,
                                  cbind(grid[uMinIndex[!duplicatePos], uname],
                                        grid[uMaxIndex[!duplicatePos], uname]))
        return(monoObjects)
    }
    genmonoAte <- function(monoObjects, type, audit = FALSE) {
        monoAte0 <- -A0[uMaxIndex, ] + A0[uMinIndex, ]
        monoAte1 <- A1[uMaxIndex, ] - A1[uMinIndex, ]
        if (is.null(dim(monoAte0))) monoAte0 <- matrix(monoAte0, nrow = 1)
        if (is.null(dim(monoAte1))) monoAte1 <- matrix(monoAte1, nrow = 1)
        duplicatePos <- duplicated(cbind(monoAte0, monoAte1))
        monoAte0 <- matrix(monoAte0[!duplicatePos, ], nrow = sum(!duplicatePos))
        monoAte1 <- matrix(monoAte1[!duplicatePos, ], nrow = sum(!duplicatePos))
        if (!audit) {
            monoAte <- cbind(matrix(0, nrow = nrow(monoAte0), ncol = 2 * sn),
                             monoAte0,
                             monoAte1)
            colnames(monoAte) <- namesA
        } else {
            monoAte <- cbind(monoAte0, monoAte1)
        }
        monoObjects$monoAte <- rbind(monoObjects$monoAte, monoAte)
        monoObjects$monotez <- c(monoObjects$monotez,
                                 replicate(nrow(monoAte), 0))
        if (type == 1) {
            monoObjects$monotes <- c(monoObjects$monotes,
                                     replicate(nrow(monoAte), ">="))
        }
        if (type == -1) {
            monoObjects$monotes <- c(monoObjects$monotes,
                                     replicate(nrow(monoAte), "<="))
        }
        monoObjects$monoAteseq <- rbind(monoObjects$monoAteseq,
                                        cbind(seq(1, nrow(monoAte)) +
                                              monoObjects$countseq,
                                              type))
        monoObjects$countseq <- monoObjects$countseq + nrow(monoAte)
        monoObjects$monomap <- c(monoObjects$monomap,
                                 gridmap[uMinIndex[!duplicatePos]])
        monoObjects$umap <- rbind(monoObjects$umap,
                                  cbind(grid[uMinIndex[!duplicatePos], uname],
                                        grid[uMaxIndex[!duplicatePos], uname]))
        return(monoObjects)
    }
    ## Implement functions---monotonicity matrices formed immediately
    ## in order to save memory. If solutions are passed, the audit is
    ## also performed, i.e. the function checks whether or not the
    ## monotonicity constraints are satisfied.
    if (!audit) {
        monoA <- NULL
    } else {
        monoA <- list()
    }
    m0.type <- 0
    m1.type <- 0
    mte.type <- 0
    if (audit) {
        monoA0IncSeq <- NULL
        monoA0DecSeq <- NULL
        monoA1IncSeq <- NULL
        monoA1DecSeq <- NULL
        monoAteIncSeq <- NULL
        monoAteDecSeq <- NULL
        diff <- NULL
        map <- NULL
        umap <- NULL
    }
    if (try(m0.inc, silent = TRUE) == TRUE) {
        monoList <- genmonoA0(monoList, 1, audit)
        m0.type <- m0.type + 1
    }
    if (try(m0.dec, silent = TRUE) == TRUE) {
        monoList <- genmonoA0(monoList, -1, audit)
        m0.type <- m0.type + 2
    }
    if (!is.null(monoList$monoA0)) {
        if (!audit) {
            monoA <- rbind(monoA, monoList$monoA0)
            rm(m0.type)
            monoList$monoA0 <- NULL
        } else {
            violateDiff <- cbind(monoList$monoA0 %*% solution.m0.min,
                                 monoList$monoA0 %*% solution.m0.max)
            violateDiff <- apply(violateDiff, 1, max)
            negatepos <- which(monoList$mono0s == ">=")
            violateDiff[negatepos] <- -violateDiff[negatepos]
            violatePos <- violateDiff > audit.tol
            if (sum(violatePos) > 0) {
                diff <- c(diff, violateDiff[violatePos])
                monoA$m0 <- matrix(monoList$monoA0[violatePos, ],
                                      nrow = sum(violatePos))
                map <- c(map, monoList$monomap[violatePos])
                umap <- c(umap, monoList$umap[violatePos, 2])
                if (m0.type == 1) {
                    monoA0IncSeq <- seq(sum(violatePos))
                    monoA0DecSeq <- NULL
                } else if (m0.type == 2) {
                    monoA0IncSeq <- NULL
                    monoA0DecSeq <- seq(sum(violatePos))
                } else if (m0.type == 3) {
                    monoA0IncSeq <- seq(sum(violatePos[violatePos <= (
                        nrow(monoList$monoA0) / 2)]))
                    monoA0DecSeq <- seq(sum(violatePos[violatePos > (
                        nrow(monoList$monoA0) / 2)]))
                }
            }
            rm(m0.type)
            monoList$monoA0 <- NULL
            monoList$monomap <- NULL
            monoList$umap <- NULL
        }
    }
    if (try(m1.inc, silent = TRUE) == TRUE) {
        monoList <- genmonoA1(monoList, 1, audit)
        m1.type <- m1.type + 1
    }
    if (try(m1.dec, silent = TRUE) == TRUE) {
        monoList <- genmonoA1(monoList, -1, audit)
        m1.type <- m1.type + 2
    }
    if (!is.null(monoList$monoA1)) {
        if (!audit) {
            monoA <- rbind(monoA, monoList$monoA1)
            monoList$monoA1 <- NULL
        } else {
            violateDiff <- cbind(monoList$monoA1 %*% solution.m1.min,
                                 monoList$monoA1 %*% solution.m1.max)
            violateDiff <- apply(violateDiff, 1, max)
            negatepos <- which(monoList$mono1s == ">=")
            violateDiff[negatepos] <- -violateDiff[negatepos]
            violatePos <- violateDiff > audit.tol
            if (sum(violatePos) > 0) {
                diff <- c(diff, violateDiff[violatePos])
                monoA$m1 <- matrix(monoList$monoA1[violatePos, ],
                                      nrow = sum(violatePos))
                map <- c(map, monoList$monomap[violatePos])
                umap <- c(umap, monoList$umap[violatePos, 2])
                if (m1.type == 1) {
                    monoA1IncSeq <- seq(sum(violatePos))
                    monoA1DecSeq <- NULL
                } else if (m1.type == 2) {
                    monoA1IncSeq <- NULL
                    monoA1DecSeq <- seq(sum(violatePos))
                } else if (m1.type == 3) {
                    monoA1IncSeq <- seq(sum(violatePos[violatePos <= (
                        nrow(monoList$monoA1) / 2)]))
                    monoA1DecSeq <- seq(sum(violatePos[violatePos > (
                        nrow(monoList$monoA1) / 2)]))
                }
            }
            rm(m1.type)
            monoList$monoA0 <- NULL
            monoList$monomap <- NULL
            monoList$umap <- NULL
        }
    }
    if (try(mte.inc, silent = TRUE) == TRUE) {
        monoList <- genmonoAte(monoList, 1, audit)
        mte.type <- mte.type + 1
    }
    if (try(mte.dec, silent = TRUE) == TRUE) {
        monoList <- genmonoAte(monoList, -1, audit)
        mte.type <- mte.type + 2
    }
    if (!is.null(monoList$monoAte)) {
        if (!audit) {
            monoA <- rbind(monoA, monoList$monoAte)
            monoList$monoAte <- NULL
        } else {
            violateDiff <- cbind(monoList$monoAte %*%
                                 c(solution.m0.min, solution.m1.min),
                                 monoList$monoAte %*%
                                 c(solution.m0.max, solution.m1.max))
            violateDiff <- apply(violateDiff, 1, max)
            negatepos <- which(monoList$monotes == ">=")
            violateDiff[negatepos] <- -violateDiff[negatepos]
            violatePos <- violateDiff > audit.tol
            if (sum(violatePos) > 0) {
                diff <- c(diff, violateDiff[violatePos])
                monoA$mte <- matrix(monoList$monoAte[violatePos, ],
                                       nrow = sum(violatePos))
                map <- c(map, monoList$monomap[violatePos])
                umap <- c(umap, monoList$umap[violatePos, 2])
                if (mte.type == 1) {
                    monoAteIncSeq <- seq(sum(violatePos))
                    monoAteDecSeq <- NULL
                } else if (mte.type == 2) {
                    monoAteIncSeq <- NULL
                    monoAteDecSeq <- seq(sum(violatePos))
                } else if (mte.type == 3) {
                    monoAteIncSeq <- seq(sum(violatePos[violatePos <= (
                        nrow(monoList$monoAte) / 2)]))
                    monoAteDecSeq <- seq(sum(violatePos[violatePos > (
                        nrow(monoList$monoAte) / 2)]))
                }
            }
            rm(mte.type)
            monoList$monoAte <- NULL
            monoList$monomap <- NULL
            monoList$umap <- NULL
        }
    }
    if (!audit) {
        ## Combine remaining vectors and return
        monos   <- c(monoList$mono0s, monoList$mono1s, monoList$monotes)
        monorhs <- c(monoList$mono0z, monoList$mono1z, monoList$monotez)
        if (!is.null(monoList$monomap)) {
            monomap <- matrix(monoList$monomap, ncol = 1)
            colnames(monomap) <- c("grid.X.index")
        } else {
            monomap <- NULL
        }
        if (!is.null(monoList$umap)) {
            umap <- monoList$umap
            colnames(umap) <- c("u1", "u2")
        } else {
            umap <- NULL
        }
        if ((hasArg(m0.inc) && m0.inc == TRUE) |
            (hasArg(m0.dec) && m0.dec == TRUE)) {
            mono0seq <- monoList$monoA0seq
            if (is.null(dim(mono0seq))) mono0seq <- matrix(mono0seq, nrow = 1)
            colnames(mono0seq) <- c("row", "type (inc+/dec-)")
        } else {
            mono0seq <- NULL
        }
        if ((hasArg(m1.inc) && m1.inc == TRUE) |
            (hasArg(m1.dec) && m1.dec == TRUE)) {
            mono1seq <- monoList$monoA1seq
            if (is.null(dim(mono1seq))) mono1seq <- matrix(mono1seq, nrow = 1)
            colnames(mono1seq) <- c("row", "type (inc+/dec-)")
        } else {
            mono1seq <- NULL
        }
        if ((hasArg(mte.inc) && mte.inc == TRUE) |
            (hasArg(mte.dec) && mte.dec == TRUE)) {
            monoteseq <- monoList$monoAteseq
            if (is.null(dim(monoteseq))) monoteseq <- matrix(monoteseq,
                                                             nrow = 1)
            colnames(monoteseq) <- c("row", "type (inc+/dec-)")
        } else {
            monoteseq <- NULL
        }
        return(list(A = monoA,
                    sense = monos,
                    rhs = monorhs,
                    map = monomap,
                    umap = umap,
                    mono0seq = mono0seq,
                    mono1seq = mono1seq,
                    monoteseq = monoteseq))
    } else {
        ## Construct violation matrix
        if (sum(unlist(lapply(monoA, function(x) nrow(x)))) > 0) {
            violateMat <- data.frame(pos = c(monoA0IncSeq,
                                             monoA0DecSeq,
                                             monoA1IncSeq,
                                             monoA1DecSeq,
                                             monoAteIncSeq,
                                             monoAteDecSeq),
                                     type = c(rep(7, length(monoA0IncSeq)),
                                              rep(8, length(monoA0DecSeq)),
                                              rep(9, length(monoA1IncSeq)),
                                              rep(10, length(monoA1DecSeq)),
                                              rep(11, length(monoAteIncSeq)),
                                              rep(12, length(monoAteDecSeq))),
                                     grid.x = map,
                                     grid.u = umap,
                                     diff = diff)
            violateMat$group.name <- paste0(violateMat$type, ".",
                                            violateMat$grid.x)
            violateMat$type.string <- factor(violateMat$type,
                                             levels = seq(7, 12),
                                             labels = c('m0.inc',
                                                        'm0.dec',
                                                        'm1.inc',
                                                        'm1.dec',
                                                        'mte.inc',
                                                        'mte.dec'))
            return(list(monoA = monoA,
                        violateMat = violateMat))
        } else {
            return(NULL)
        }
    }
}


#' Combining the boundedness and monotonicity constraint objects
#'
#' This function simply combines the objects associated with the
#' boundedness constraints and the monotonicity constraints.
#' @param bdA list containing the constraint matrix, vector of
#'     inequalities, and RHS vector associated with the boundedness
#'     constraints.
#' @param monoA list containing the constraint matrix, vector on
#'     inequalities, and RHS vector associated with the monotonicity
#'     constraints.
#' @return a list containing a unified constraint matrix, unified
#'     vector of inequalities, and unified RHS vector for the
#'     boundedness and monotonicity constraints of an LP/QCQP problem.
combinemonobound <- function(bdA, monoA) {
    mbA    <- NULL
    mbs    <- NULL
    mbrhs  <- NULL
    mbmap  <- NULL
    mbumap <- NULL
    if (!is.null(bdA)) {
        mbA    <- rbind(mbA, bdA$A)
        mbs    <- c(mbs, bdA$sense)
        mbrhs  <- c(mbrhs, bdA$rhs)
        mbmap  <- c(mbmap, bdA$map)
        mbumap <- rbind(mbumap, cbind(bdA$umap, bdA$umap))
        ## bdA$umap is cbind'ed twice to be conformable with
        ## monoA$umap, where we must keep track of pairs of u's.
    }
    if (!is.null(monoA)) {
        mbA      <- rbind(mbA, monoA$A)
        mbs      <- c(mbs, monoA$sense)
        mbrhs    <- c(mbrhs, monoA$rhs)
        mbmap    <- c(mbmap, monoA$map)
        mbumap   <- rbind(mbumap, monoA$umap)
    }
    return(list(mbA = mbA,
                mbs = mbs,
                mbrhs  = mbrhs,
                mbmap  = mbmap,
                mbumap = mbumap))
}

#' Generating monotonicity and boundedness constraints
#'
#' This is a wrapper function generating the matrices and vectors
#' associated with the monotonicity and boundedness constraints
#' declared by the user. Since this function generates all the
#' components required for the shape constraints, it is also the
#' function that performs the audit. That is, MTR coefficients are
#' passed, then this function will verify whether they satisfy the
#' shape constraints.
#' @param pm0 A list of the monomials in the MTR for d = 0.
#' @param pm1 A list of the monomials in the MTR for d = 1.
#' @param support a matrix for the support of all variables that enter
#'     into the MTRs.
#' @param grid_index a vector, the row numbers of \code{support} used
#'     to generate the grid preceding the audit.
#' @param uvec a vector, the points in the interval [0, 1] that the
#'     unobservable takes on.
#' @param splinesobj a list of lists. Each of the inner lists contains
#'     details on the splines declared in the MTRs.
#' @param monov name of variable for which the monotonicity conditions
#'     applies to.
#' @param uname name declared by user to represent the unobservable
#'     term in the MTRs.
#' @param m0 one-sided formula for marginal treatment response
#'     function for the control group. The formula may differ from
#'     what the user originally input in \code{\link{ivmte}}, as the
#'     spline components should have been removed. This formula is
#'     simply a linear combination of all covariates that enter into
#'     the original \code{m0} declared by the user in
#'     \code{\link{ivmte}}.
#' @param m1 one-sided formula for marginal treatment response
#'     function for the treated group. The formula may differ from
#'     what the user originally input in \code{\link{ivmte}}, as the
#'     spline components should have been removed. This formula is
#'     simply a linear combination of all covariates that enter into
#'     the original \code{m1} declared by the user in
#'     \code{\link{ivmte}}.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set.
#' @param gstar0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param gstar1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param m0.lb scalar, lower bound on MTR for control group.
#' @param m0.ub scalar, upper bound on MTR for control group.
#' @param m1.lb scalar, lower bound on MTR for treated group.
#' @param m1.ub scalar, upper bound on MTR for treated group.
#' @param mte.lb scalar, lower bound on MTE.
#' @param mte.ub scalar, upper bound on MTE.
#' @param m0.dec boolean, indicating whether the MTR for the control
#'     group is monotone decreasing.
#' @param m0.inc boolean, indicating whether the MTR for the control
#'     group is monotone increasing.
#' @param m1.dec boolean, indicating whether the MTR for the treated
#'     group is monotone decreasing.
#' @param m1.inc boolean, indicating whether the MTR for the treated
#'     group is monotone increasing.
#' @param mte.dec boolean, indicating whether the MTE is monotone
#'     decreasing.
#' @param mte.inc boolean, indicating whether the MTE is monotone
#'     increasing.
#' @param solution.m0.min vector, the coefficients for the MTR for
#'     \code{D = 0} corresponding to the lower bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m1.min vector, the coefficients for the MTR for
#'     \code{D = 1} corresponding to the lower bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m0.max vector, the coefficients for the MTR for
#'     \code{D = 0} corresponding to the upper bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param solution.m1.max vector, the coefficients for the MTR for
#'     \code{D = 1} corresponding to the upper bound of the target
#'     parameter. If passed, this will initiate checks of shape
#'     constraints.
#' @param audit.tol feasibility tolerance when performing the
#'     audit. By default to set to be equal \code{1e-06}. This
#'     parameter should only be changed if the feasibility tolerance
#'     of the solver is changed, or if numerical issues result in
#'     discrepancies between the solver's feasibility check and the
#'     audit.
#' @param qp boolean, set to \code{TRUE} if the direct MTR
#'     regression via QP is used.
#' @return a list containing a unified constraint matrix, unified
#'     vector of inequalities, and unified RHS vector for the
#'     boundedness and monotonicity constraints of an LP/QCQP problem.
genmonoboundA <- function(pm0, pm1, support, grid_index, uvec,
                          splinesobj, monov, uname, m0, m1, sset,
                          gstar0, gstar1, m0.lb, m0.ub, m1.lb, m1.ub,
                          mte.lb, mte.ub, m0.dec, m0.inc, m1.dec,
                          m1.inc, mte.dec, mte.inc,
                          solution.m0.min = NULL,
                          solution.m1.min = NULL,
                          solution.m0.max = NULL,
                          solution.m1.max = NULL, audit.tol,
                          qp) {
    if (!is.null(solution.m0.min) && !is.null(solution.m1.min) &&
        !is.null(solution.m0.max) && !is.null(solution.m1.max)) {
        audit <- TRUE
    } else {
        audit <- FALSE
    }
    call <- match.call()
    splines <- list(splinesobj[[1]]$splineslist,
                    splinesobj[[2]]$splineslist)
    splinesinter <- list(splinesobj[[1]]$splinesinter,
                         splinesobj[[2]]$splinesinter)
    if (length(grid_index) == 0) {
        noX <- TRUE
    } else {
        noX <- FALSE
    }
    ## First construct the non-spline U terms (to check for redundant
    ## points in the U grid)
    u0mat <- NULL
    if (!is.null(pm0)) {
        uExp <- unique(sort(pm0$exporder))
        uExp <- uExp[uExp > 0]
        if (length(uExp) > 0) {
            for (i in uExp) {
                u0mat <- cbind(u0mat, uvec ^ i)
            }
            colnames(u0mat) <- paste0("u0.", uExp)
        }
    }
    u1mat <- NULL
    if (!is.null(pm1)) {
        uExp <- unique(sort(pm1$exporder))
        uExp <- uExp[uExp > 0]
        if (length(uExp) > 0) {
            for (i in uExp) {
                u1mat <- cbind(u1mat, uvec ^ i)
            }
            colnames(u1mat) <- paste0("u1.", uExp)
        }
    }
    ## Now construct the splines U terms (to check for redundant
    ## points in the U grid)
    us0mat <- Reduce(cbind,
                     genBasisSplines(splines = splines[[1]],
                                     x = uvec,
                                     d = 0))
    us1mat <- Reduce(cbind,
                     genBasisSplines(splines = splines[[2]],
                                     x = uvec,
                                     d = 1))
    fullU0mat <- cbind(u0mat, us0mat)
    fullU1mat <- cbind(u1mat, us1mat)
    u0unique <- uvec[!duplicated(fullU0mat)]
    u1unique <- uvec[!duplicated(fullU1mat)]
    ## Now update U grid to only include non-redundant values
    uvec <- sort(unique(c(u0unique, u1unique)))
    if (length(uvec) == 0) uvec <- 0
    rm(u0mat, u1mat,
       us0mat, us1mat,
       fullU0mat, fullU1mat,
       u0unique, u1unique)
    ## Generate the first iteration of the grid
    if (noX) {
        grid <- data.frame(uvec)
        grid$.grid.order <- seq(1, length(uvec))
        grid$.u.order <- seq(1, length(uvec))
        colnames(grid) <- c(uname, ".grid.order", ".u.order")
        grid_index <- rownames(grid)
        gridobj <- list(grid = grid,
                        map  = replicate(length(uvec), 0))
    } else {
        gridobj <- gengrid(grid_index,
                           support,
                           uvec,
                           uname)
    }
    if (is.null(splines[[1]]) && is.null(splines[[2]])) {
        A0 <- design(formula = m0, data = gridobj$grid)$X
        A1 <- design(formula = m1, data = gridobj$grid)$X
        A0 <- cbind(A0,
                    .grid.order = gridobj$grid$.grid.order,
                    .u.order = gridobj$grid$.u.order)
        A1 <- cbind(A1,
                    .grid.order = gridobj$grid$.grid.order,
                    .u.order = gridobj$grid$.u.order)
    } else {
        ## Construct base A0 and A1 matrices using the non-spline
        ## formulas. The variables interacting with the splines will
        ## be appended to this.
        if (is.null(m0)) {
            A0 <- design(formula = as.formula(paste("~ 0 +", uname)),
                         data = gridobj$grid)$X
        } else {
            m0 <- as.formula(paste(gsub("\\s+", " ",
                                        Reduce(paste, deparse(m0))), "+",
                                   uname))
            A0 <- design(formula = m0, data = gridobj$grid)$X
            colnames(A0) <- parenthBoolean(colnames(A0))
        }
        if (is.null(m1)) {
            A1 <- design(formula = as.formula(paste("~ 0 +", uname)),
                         data = gridobj$grid)$X
        } else {
            m1 <- as.formula(paste(gsub("\\s+", " ",
                                        Reduce(paste, deparse(m1))), "+",
                                   uname))
            A1 <- design(formula = m1, data = gridobj$grid)$X
            colnames(A1) <- parenthBoolean(colnames(A1))
        }
        dlist <- NULL
        if (!is.null(splines[[1]])) dlist <- c(dlist, 0)
        if (!is.null(splines[[2]])) dlist <- c(dlist, 1)
        for (d in dlist) {
            nonSplinesDmat <- NULL
            splinesD <- splines[[d + 1]]
            for (j in 1:length(splinesD)) {
                for (k in 1:length(splinesD[[j]])) {
                    if (splinesD[[j]][k] != "1") {
                        ## Check if variable is a factor and whether
                        ## it has sufficient variance
                        if (substr(splinesD[[j]][k], 1, 7) == "factor(") {
                            tmpVar <- substr(splinesD[[j]][k],
                                             8,
                                             nchar(splinesD[[j]][k]) - 1)
                            tmpVarSupport <- unique(gridobj$grid[, tmpVar])
                            if (length(tmpVarSupport) == 1) {
                                tmpDmat <-
                                    matrix(rep(1, nrow(gridobj$grid)), ncol = 1)
                                colnames(tmpDmat) <- paste0("factor(",
                                                            tmpVar,
                                                            ")", tmpVarSupport)
                            }
                        }
                        if (!exists("tmpDmat")) {
                            tmpDmat <-
                                design(as.formula(paste("~ 0 +",
                                                        splinesD[[j]][k])),
                                       gridobj$grid)$X
                        }
                        nonSplinesDmat <- cbind(nonSplinesDmat, tmpDmat)
                        rm(tmpDmat)
                    } else {
                        nonSplinesDmat <- cbind(nonSplinesDmat,
                                                design(~ 1, gridobj$grid)$X)
                    }
                }
                colnames(nonSplinesDmat) <-
                    parenthBoolean(colnames(nonSplinesDmat))
            }
            ## Only keep the variables that are not already in A0 or A1
            currentNames <- colnames(nonSplinesDmat)
            keepPos <- ! colnames(nonSplinesDmat) %in%
                colnames(get(paste0("A", d)))
            nonSplinesDmat <- as.matrix(nonSplinesDmat[, keepPos])
            colnames(nonSplinesDmat) <- currentNames[keepPos]
            assign(paste0("A", d),
                   cbind(get(paste0("A", d)), nonSplinesDmat))
        }
        A0 <- cbind(A0,
                    .grid.order = gridobj$grid$.grid.order,
                    .u.order = gridobj$grid$.u.order)
        A1 <- cbind(A1,
                    .grid.order = gridobj$grid$.grid.order,
                    .u.order = gridobj$grid$.u.order)
        basisList <- list(genBasisSplines(splines = splines[[1]],
                                          x = uvec,
                                          d = 0),
                          genBasisSplines(splines = splines[[2]],
                                          x = uvec,
                                          d = 1))

        ## Generate interaction with the splines.
        ## Indexing in the loops takes the following structure:
        ## j: splines index
        ## v: interaction index
        colnames(A0) <- parenthBoolean(colnames(A0))
        colnames(A1) <- parenthBoolean(colnames(A1))
        namesA0 <- colnames(A0)
        namesA1 <- colnames(A1)
        namesA0length <- sapply(namesA0, function(x) {
            length(unlist(strsplit(x, ":")))
        })
        namesA1length <- sapply(namesA1, function(x) {
            length(unlist(strsplit(x, ":")))
        })
        for (d in 0:1) {
            namesA <- get(paste0("namesA", d))
            namesAlength <- get(paste0("namesA", d, "length"))
            if (!is.null(basisList[[d + 1]])) {
                for (j in 1:length(splines[[d + 1]])) {
                    for (v in 1:length(splines[[d + 1]][[j]])) {
                        bmat <- cbind(uvec, basisList[[d + 1]][[j]])
                        colnames(bmat)[1] <- uname
                        iName <- splines[[d + 1]][[j]][v]
                        if (iName != "1") {
                            isFactor  <- FALSE
                            isBoolean <- FALSE
                            iNamePos  <- NULL
                            iNameList <- unlist(strsplit(iName, ":"))
                            namesAscore <- rep(0, times = length(namesA))
                            for (q in iNameList) {
                                if (substr(q, nchar(q), nchar(q)) == ")") {
                                    q <- gsub("\\)", "\\\\)",
                                              gsub("\\(", "\\\\(", q))
                                    q <- gsub("\\+", "\\\\+", q)
                                    q <- gsub("\\*", "\\\\*", q)
                                    q <- gsub("\\^", "\\\\^", q)
                                    namesApos <-
                                        as.integer(
                                            grepl(paste0(q, "[0-9A-Za-z._]*"),
                                                  namesA))
                                } else {
                                    namesApos <- as.integer(namesA == q)
                                }
                                namesAscore <- namesAscore + namesApos
                            }
                            iNamePos <- (namesAscore == length(iNameList)) *
                                (namesAscore == namesAlength)
                            iNamePos <- which(iNamePos == 1)
                            for (r in iNamePos) {
                                bmatTmp <-
                                    merge(
                                        get(paste0("A", d))[, c(uname,
                                                                namesA[r],
                                                                ".grid.order")],
                                        bmat, by = uname)
                                tmpX <- bmatTmp[, 4:ncol(bmatTmp)]
                                if (is.null(dim(tmpX))) {
                                    tmpX <- matrix(tmpX, ncol = 1)
                                }
                                bmatTmp[, 4:ncol(bmatTmp)] <-
                                    sweep(x = tmpX,
                                          MARGIN = 1,
                                          STATS = bmatTmp[, namesA[r]],
                                          FUN = "*")
                                namesB <-
                                    paste0(colnames(bmatTmp)[4:ncol(bmatTmp)],
                                           ":", namesA[r])
                                colnames(bmatTmp)[4:ncol(bmatTmp)] <- namesB
                                newA <- merge(get(paste0("A", d)),
                                              bmatTmp[, c(".grid.order",
                                                          namesB)],
                                              by = ".grid.order")
                                assign(paste0("A", d), newA)
                                rm(bmatTmp, newA, tmpX)
                            }
                        } else {
                            namesA <- colnames(get(paste0("A", d)))
                            namesB <- paste0(colnames(bmat)[2:ncol(bmat)],
                                             ":", iName)
                            colnames(bmat)[2:ncol(bmat)] <- namesB
                            newA <- merge(get(paste0("A", d)),
                                          bmat,
                                          by = uname)
                            assign(paste0("A", d), newA)
                            rm(newA)
                        }
                        rm(bmat)
                    }
                }
            }
        }
        rownames(A0) <- A0[, ".grid.order"]
        rownames(A1) <- A1[, ".grid.order"]
    }
    A0 <- A0[order(A0[, ".grid.order"]), ]
    A1 <- A1[order(A1[, ".grid.order"]), ]
    ## If both m0 and m1 are just the intercepts, then the objects A0
    ## and A1 will not be matrices. This will cause an error when
    ## trying to account for column names. So convert A0 and A1 into
    ## matrices.
    if (is.null(dim(A0))) {
        tmpNames <- names(A0)
        A0 <- matrix(A0, nrow = 1)
        colnames(A0) <- tmpNames
        rownames(A0) <- 1
        rm(tmpNames)
    }
    if (is.null(dim(A1))) {
        tmpNames <- names(A1)
        A1 <- matrix(A1, nrow = 1)
        colnames(A1) <- tmpNames
        rownames(A1) <- 1
        rm(tmpNames)
    }
    ## Rename columns so they match with the names in vectors gstar0
    ## and gstar1 (the problem stems from the unpredictable ordering
    ## of variables in interaction terms).
    for (d in c(0, 1)) {
        Amat <- get(paste0("A", d))
        if (d == 0) rm(A0)
        if (d == 1) rm(A1)
        gvec <- get(paste0("gstar", d))
        Apos <- NULL
        failTerms <- which(!names(gvec) %in% colnames(Amat))
        for (fail in failTerms) {
            vars <- strsplit(names(gvec)[fail], ":")[[1]]
            varsPerm <- permute(vars)
            varsPerm <- unlist(lapply(varsPerm,
                                      function(x) paste(x, collapse = ":")))
            correctPos <- unique(which(varsPerm %in% colnames(Amat)))
            if (length(correctPos) > 0) {
                Aname <- varsPerm[correctPos]
                Apos <- c(Apos, which(colnames(Amat) == Aname))
            }
        }
        colnames(Amat)[Apos] <- names(gvec)[failTerms]
        assign(paste0("A", d), Amat)
        rm(Amat)
    }
    ## Some columns maybe missing relative to gstar0/gstar1 becuase
    ## the grid is not large enough, and so does not contain all
    ## factor variables. Fill these columns with 0
    missingA0 <- !(names(gstar0) %in% colnames(A0))
    A0 <- as.data.frame(A0)
    for (i in names(gstar0)[missingA0]) {
        A0[, i] <- 0
    }
    missingA1 <- !(names(gstar1) %in% colnames(A1))
    A1 <- as.data.frame(A1)
    for (i in names(gstar1)[missingA1]) {
        A1[, i] <- 0
    }
    ## keep only the columns that are in the MTRs (A0 and A1 matrices
    ## potentially include extraneous columns)
    A0 <- as.matrix(A0[, names(gstar0)])
    A1 <- as.matrix(A1[, names(gstar1)])
    colnames(A0) <- names(gstar0)
    colnames(A1) <- names(gstar1)
    ## Remove duplicate rows---these become redundant in the audit
    duplicatePos <- duplicated(cbind(A0, A1))
    namesA0 <- colnames(A0)
    namesA1 <- colnames(A1)
    A0 <- matrix(A0[!duplicatePos, ], nrow = sum(!duplicatePos))
    A1 <- matrix(A1[!duplicatePos, ], nrow = sum(!duplicatePos))
    colnames(A0) <- namesA0
    colnames(A1) <- namesA1
    rm(namesA0, namesA1)
    gridobj$grid <- gridobj$grid[!duplicatePos, ]
    gridobj$map <- gridobj$map[!duplicatePos]
    ## Construct placeholders
    bdA <- NULL
    monoA <- NULL
    lb0seq <- NULL
    lb1seq <- NULL
    lbteseq <- NULL
    ub0seq <- NULL
    ub1seq <- NULL
    ubteseq <- NULL
    mono0seq <- NULL
    mono1seq <- NULL
    monomteseq <- NULL
    ## generate matrices for imposing bounds on m0 and m1 and
    ## treatment effects
    if (hasArg(m0.lb) | hasArg(m0.ub) |
        hasArg(m1.lb) | hasArg(m1.lb) |
        hasArg(mte.lb) | hasArg(mte.ub)) {
        boundlist  <- c("uname",
                        "m0.lb", "m0.ub",
                        "m1.lb", "m1.ub",
                        "mte.lb", "mte.ub",
                        "solution.m0.min", "solution.m1.min",
                        "solution.m0.max", "solution.m1.max",
                        "audit.tol", "qp")
        boundAcall <- modcall(call,
                              newcall = genboundA,
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
        if (!(length(uvec) == 1 && uvec == 0)) {
            monolist  <- c("m0.dec", "m0.inc",
                           "m1.dec", "m1.inc",
                           "mte.dec", "mte.inc",
                           "solution.m0.min", "solution.m1.min",
                           "solution.m0.max", "solution.m1.max",
                           "audit.tol", "qp")
            monoAcall <- modcall(call,
                                 newcall = genmonoA,
                                 keepargs = monolist,
                                 newargs = list(A0 = quote(A0),
                                                A1 = quote(A1),
                                                sset = quote(sset),
                                                gridobj = quote(gridobj),
                                                uname = uname,
                                                gstar0 = quote(gstar0),
                                                gstar1 = quote(gstar1)))
            monoA <- eval(monoAcall)
        }
    }
    if (audit) {
        ## Return violation matrices, if the function is used to
        ## perform the audit
        return(list(bounds = bdA,
                    mono = monoA))
    }
    ## Update bound sequence counts
    if (!is.null(bdA$lb0seq)) lb0seq <- bdA$lb0seq
    if (!is.null(bdA$lb1seq)) lb1seq <- bdA$lb1seq
    if (!is.null(bdA$lbteseq)) lbteseq <- bdA$lbteseq
    if (!is.null(bdA$ub0seq)) ub0seq <- bdA$ub0seq
    if (!is.null(bdA$ub1seq)) ub1seq <- bdA$ub1seq
    if (!is.null(bdA$ubteseq)) ubteseq <- bdA$ubteseq
    ## Update monotonicity sequence counts
    boundLength <- length(lb0seq) + length(lb1seq) + length(lbteseq) +
        length(ub0seq) + length(ub1seq) + length(ubteseq)
    if (!is.null(monoA$mono0seq)) {
        mono0seq <- monoA$mono0seq
        monoA$mono0seq <- NULL
        mono0seq[, 1] <- mono0seq[, 1] + boundLength
    }
    if (!is.null(monoA$mono1seq)) {
        mono1seq <- monoA$mono1seq
        monoA$mono1seq <- NULL
        mono1seq[, 1] <- mono1seq[, 1] + boundLength
    }
    if (!is.null(monoA$monoteseq)) {
        monoteseq <- monoA$monoteseq
        monoA$monoteseq <- NULL
        monoteseq[, 1] <- monoteseq[, 1] + boundLength
    }
    mbA <- list(bdA$A, monoA$A)
    mbs <- list(bdA$sense, monoA$sense)
    mbrhs <- list(bdA$rhs, monoA$rhs)
    mbmap <- list(bdA$map, monoA$map)
    mbumap <- list(bdA$umap, ## This needs to be doubled/cbinded
                   monoA$umap)
    rm(bdA, monoA)
    mbA <- Reduce(rbind, mbA)
    mbs <- Reduce(c, mbs)
    mbrhs <- Reduce(c, mbrhs)
    mbmap <- Reduce(c, mbmap)
    mbumap <- rbind(cbind(mbumap[[1]], mbumap[[1]]),
                    mbumap[[2]])
    output <- list(mbA = mbA,
                   mbs = mbs,
                   mbrhs  = mbrhs)
                   ## mbmap  = mbmap,
                   ## mbumap = mbumap)
    rm(mbA, mbs, mbrhs, mbmap, mbumap)
    ## output$gridobj <- gridobj
    output$lb0seq  <- lb0seq
    output$lb1seq  <- lb1seq
    output$lbteseq <- lbteseq
    output$ub0seq  <- ub0seq
    output$ub1seq  <- ub1seq
    output$ubteseq <- ubteseq
    ## Set rownames of constraint matrix to describe constraint
    rownames(output$mbA) <- rep('', nrow(output$mbA))
    if (length(lb0seq) > 0) rownames(output$mbA)[lb0seq] <-
                                 rep('lb0', length(lb0seq))
    if (length(lb1seq) > 0) rownames(output$mbA)[lb1seq] <-
                                 rep('lb1', length(lb1seq))
    if (length(lbteseq) > 0) rownames(output$mbA)[lbteseq] <-
                                 rep('lbte', length(lbteseq))
    if (length(ub0seq) > 0) rownames(output$mbA)[ub0seq] <-
                                 rep('ub0', length(ub0seq))
    if (length(ub1seq) > 0) rownames(output$mbA)[ub1seq] <-
                                 rep('ub1', length(ub1seq))
    if (length(ubteseq) > 0) rownames(output$mbA)[ubteseq] <-
                                 rep('ubte', length(ubteseq))
    if (exists("mono0seq")) {
        output$mono0seq <- mono0seq
        rownames(output$mbA)[mono0seq]  <- 'mono0'
        rm(mono0seq)
    }
    if (exists("mono1seq")) {
        output$mono1seq <- mono1seq
        rownames(output$mbA)[mono1seq]  <- 'mono1'
        rm(mono1seq)
    }
    if (exists("monoteseq")) {
        output$monoteseq <- monoteseq
        rownames(output$mbA)[monoteseq]  <- 'monote'
        rm(monoteseq)
    }
    return(output)
}
