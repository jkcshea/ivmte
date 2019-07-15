#' Generating the grid for the audit procedure
#'
#' This function takes in a matrix summarizing the support of the
#' covariates, as well as an evenly spaced set of points summarizing
#' the support of the unobservable variable. A Cartesian product of
#' the subset of the support of the covariates and the points in the
#' support of the unobservable generates the grid that is used for the
#' audit procedure.
#'
#' @param index a vector whose elements indicate the rows in the
#'     matrix \code{xsupport} to include in the grid.
#' @param xsupport a matrix containing all the unique combinations of
#'     the covariates included in the MTRs.
#' @param usupport a vector of evenly spaced points in the interval
#'     [0, 1], including 0 and 1. The number of points is decided by
#'     the user.
#' @param uname name declared by user to represent the unobservable
#'     term.
#' @return a list containing the grid used in the audit; a vector
#'     mapping the elements in the support of the covariates to
#'     \code{index}.
gengrid <- function(index, xsupport, usupport, uname) {

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

#' Generating the LP constraint matrix for bounds
#'
#' This function generates the component of the constraint matrix in
#' the LP problem pertaining to bounds on the MTRs and MTEs. These
#' bounds are declared by the user.
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
#' @return a constraint matrix for the LP problem, the associated
#'     vector of inequalities, and the RHS vector in the inequality
#'     constraint. The objects pertain only to the boundedness
#'     constraints declared by the user.
genboundA <- function(A0, A1, sset, gridobj, uname,
                          m0.lb, m0.ub, m1.lb, m1.ub, mte.lb, mte.ub) {

    sn <- length(sset)
    grid <- gridobj$grid
    gridmap <- gridobj$map

    namesA0 <- colnames(A0)
    namesA1 <- colnames(A1)
    namesA  <- c(seq(1, 2 * sn),
                 namesA0,
                 namesA1)

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

    ## Generate matrices for imposing bounds on m0
    if (hasArg(m0.ub) | hasArg(m0.lb)) {

        bdA0 <- cbind(matrix(0, nrow = nrow(grid), ncol = 2 * sn),
                      A0,
                      matrix(0, nrow = nrow(A1), ncol = ncol(A1)))

        colnames(bdA0) <- namesA

        if (is.numeric(try(m0.ub, silent = TRUE))) {
            ubdA0 <- bdA0
            m0ub  <- replicate(nrow(A0), m0.ub)
            m0ubs <- replicate(nrow(A0), "<=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            ubdA0seq <- seq(1, nrow(A0))
        }
        if (is.numeric(try(m0.lb, silent = TRUE))) {
            lbdA0 <- bdA0
            m0lb  <- replicate(nrow(A0), m0.lb)
            m0lbs <- replicate(nrow(A0), ">=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            lbdA0seq <- seq(1, nrow(A0))
        }
    }

    ## Generate matrices for imposing bounds on m1
    if (hasArg(m1.ub) | hasArg(m1.lb)) {
        bdA1 <- cbind(matrix(0, nrow = nrow(grid), ncol = 2 * sn),
                      matrix(0, nrow = nrow(A0),   ncol = ncol(A0)),
                      A1)
        colnames(bdA1) <- namesA

        if (is.numeric(try(m1.ub, silent = TRUE))) {
            ubdA1 <- bdA1
            m1ub  <- replicate(nrow(A1), m1.ub)
            m1ubs <- replicate(nrow(A1), "<=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            ubdA1seq <- seq(1, nrow(A1))
        }
        if (is.numeric(try(m1.lb, silent = TRUE))) {
            lbdA1 <- bdA1
            m1lb  <- replicate(nrow(A1), m1.lb)
            m1lbs <- replicate(nrow(A1), ">=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            lbdA1seq <- seq(1, nrow(A1))
        }
    }
    ## Generate matrices for imposing bounds on m1 - m0
    if(hasArg(mte.lb) | hasArg(mte.ub)) {
        bdAte <- cbind(matrix(0, nrow = nrow(grid), ncol = 2 * sn),
                       -A0, A1)
        colnames(bdAte) <- namesA

        if (is.numeric(try(mte.ub, silent = TRUE))) {
            ubdAte <- bdAte
            teub  <- replicate(nrow(A1), mte.ub)
            teubs <- replicate(nrow(A1), "<=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            ubdAteseq <- seq(1, nrow(A1))
        }
        if (is.numeric(try(mte.lb, silent = TRUE))) {
            lbdAte <- bdAte
            telb  <- replicate(nrow(A1), mte.lb)
            telbs <- replicate(nrow(A1), ">=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
            lbdAteseq <- seq(1, nrow(A1))
        }
    }

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

    ## Combine matrices and return
    bdA <- rbind(lbdA0,  lbdA1,  lbdAte,
                 ubdA0,  ubdA1,  ubdAte)
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
}

#' Generate LP components of the monotonicity constraints
#'
#' This function generates the matrix and vectors associated with the
#' monotonicity constraints declared by the user. It takes in a grid
#' of the covariates on which we define the LP constraints, and then
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
#' @return constraint matrix for the LP problem. The matrix pertains
#'     only to the monotonicity conditions on the MTR and MTE declared
#'     by the user.
genmonoA <- function(A0, A1, sset, uname, gridobj, gstar0, gstar1,
                     m0.dec, m0.inc, m1.dec, m1.inc, mte.dec,
                     mte.inc) {

    un <- length(unique(gridobj$grid[, uname]))

    ## Construct index for calculating first differences
    uMaxIndex <- seq(1, nrow(A0))[-seq(from = 1, to = nrow(A0), by = un)]
    uMinIndex <- seq(1, nrow(A0))[-seq(from = un, to = nrow(A0), by = un)]

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
    ## This matrix should include all the additions 0s on the left
    ## columns
    sn <- length(sset)
    namesA0 <- colnames(A0)
    namesA1 <- colnames(A1)
    namesA  <- c(seq(1, 2 * sn),
                 namesA0,
                 namesA1)

    ## The functions below generate the constraint matrix, the sense
    ## vector, and the RHS vector associated with the monotonicity
    ## constraints for m0, m1, and the mte. In addition, mappings to
    ## the grid index and U index are constructed; also, a list of the
    ## direction of monotonicity constriants is generated (since you
    ## allow m0, m1, and mte to face both increasing and decreasing
    ## monotonicity constraints simultaneously---forcing them to be
    ## constants).
    genmonoA0 <- function(monoObjects, type) {
        monoA0 <- A0[uMaxIndex, ] - A0[uMinIndex, ]
        if (is.null(dim(monoA0))) monoA0 <- matrix(monoA0, nrow = 1)
        monoA0 <- cbind(matrix(0, nrow = nrow(monoA0), ncol = 2 * sn),
                        monoA0,
                        matrix(0, nrow = nrow(monoA0), ncol = ncol(A1)))
        colnames(monoA0) <- namesA

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
                                 gridobj$map[uMinIndex])
        monoObjects$umap <- rbind(monoObjects$umap,
                                  cbind(gridobj$grid[uMinIndex, uname],
                                        gridobj$grid[uMaxIndex, uname]))
        return(monoObjects)
    }

    genmonoA1 <- function(monoObjects, type) {
        monoA1 <- A1[uMaxIndex, ] - A1[uMinIndex, ]
        if (is.null(dim(monoA1))) monoA1 <- matrix(monoA1, nrow = 1)
        monoA1 <- cbind(matrix(0, nrow = nrow(monoA1), ncol = 2 * sn),
                        matrix(0, nrow = nrow(monoA1), ncol = ncol(A1)),
                        monoA1)
        colnames(monoA1) <- namesA

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
                                 gridobj$map[uMinIndex])
        monoObjects$umap <- rbind(monoObjects$umap,
                                  cbind(gridobj$grid[uMinIndex, uname],
                                        gridobj$grid[uMaxIndex, uname]))
        return(monoObjects)
    }

    genmonoAte <- function(monoObjects, type) {
        monoAte0 <- -A0[uMaxIndex, ] + A0[uMinIndex, ]
        monoAte1 <- A1[uMaxIndex, ] - A1[uMinIndex, ]
        if (is.null(dim(monoAte0))) monoAte0 <- matrix(monoAte0, nrow = 1)
        if (is.null(dim(monoAte1))) monoAte1 <- matrix(monoAte1, nrow = 1)
        monoAte <- cbind(matrix(0, nrow = nrow(monoAte0), ncol = 2 * sn),
                         monoAte0,
                         monoAte1)
        colnames(monoAte) <- namesA

        monoObjects$monoAte <- rbind(monoObjects$monoAte, monoAte)
        monoObjects$monotez <- c(monoObjects$monotez, replicate(nrow(monoAte), 0))
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
                                 gridobj$map[uMinIndex])
        monoObjects$umap <- rbind(monoObjects$umap,
                                  cbind(gridobj$grid[uMinIndex, uname],
                                        gridobj$grid[uMaxIndex, uname]))
        return(monoObjects)
    }

    ## Implement functions
    if (try(m0.inc, silent = TRUE) == TRUE) {
        monoList <- genmonoA0(monoList, 1)
    }
    if (try(m0.dec, silent = TRUE) == TRUE) {
        monoList <- genmonoA0(monoList, -1)
    }
    if (try(m1.inc, silent = TRUE) == TRUE) {
        monoList <- genmonoA1(monoList, 1)
    }
    if (try(m1.dec, silent = TRUE) == TRUE) {
        monoList <- genmonoA1(monoList, -1)
    }
    if (try(mte.inc, silent = TRUE) == TRUE) {
        monoList <- genmonoAte(monoList, 1)
    }
    if (try(mte.dec, silent = TRUE) == TRUE) {
        monoList <- genmonoAte(monoList, -1)
    }

    ## Combine matrices and return
    monoA <- rbind(monoList$monoA0, monoList$monoA1, monoList$monoAte)
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
        mono0seq = monoList$monoA0seq
        if (is.null(dim(mono0seq))) mono0seq <- matrix(mono0seq, nrow = 1)
        colnames(mono0seq) <- c("row", "type (inc+/dec-)")
    } else {
        mono0seq <- NULL
    }
    if ((hasArg(m1.inc) && m1.inc == TRUE) |
        (hasArg(m1.dec) && m1.dec == TRUE)) {
        mono1seq = monoList$monoA1seq
        if (is.null(dim(mono1seq))) mono1seq <- matrix(mono1seq, nrow = 1)
        colnames(mono1seq) <- c("row", "type (inc+/dec-)")
    } else {
        mono1seq <- NULL
    }

    if ((hasArg(mte.inc) && mte.inc == TRUE) |
        (hasArg(mte.dec) && mte.dec == TRUE)) {
        monoteseq = monoList$monoAteseq
        if (is.null(dim(monoteseq))) monoteseq <- matrix(monoteseq, nrow = 1)
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
#'     boundedness and monotonicity constraints of an LP problem.
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
#' declared by the user.
#' @param support a matrix for the support of all variables that enter
#'     into the MTRs.
#' @param grid_index a vector, the row numbers of \code{support} used
#'     to generate the grid preceding the audit.
#' @param uvec a vector, the points in the interval [0, 1] that the
#'     unobservable takes on.
#' @param splines a list of lists. Each of the inner lists contains
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
#' @return a list containing a unified constraint matrix, unified
#'     vector of inequalities, and unified RHS vector for the
#'     boundedness and monotonicity constraints of an LP problem.
genmonoboundA <- function(support, grid_index, uvec, splines, monov,
                          uname, m0, m1, sset, gstar0, gstar1,
                          m0.lb, m0.ub, m1.lb, m1.ub, mte.lb, mte.ub,
                          m0.dec, m0.inc, m1.dec, m1.inc, mte.dec, mte.inc) {
    call <- match.call()
    if (is.null(grid_index)) {
        noX <- TRUE
    } else {
        noX <- FALSE
    }

    ## generate the first iteration of the grid
    if (noX) {
        grid <- data.frame(uvec)
        grid$.grid.order <- seq(1, length(uvec))
        grid$.u.order <- seq(1, length(uvec))
        colnames(grid) <- c(uname, ".grid.order", ".u.order")
        grid_index <- rownames(grid)
        gridobj <- list(grid = grid,
                        map  = replicate(length(uvec), 1))
    } else {
        gridobj <- gengrid(grid_index,
                           support,
                           uvec,
                           uname)
    }
    if (is.null(splines[[1]]) & is.null(splines[[2]])) {
        A0 <- design(formula = m0, data = gridobj$grid)$X
        A1 <- design(formula = m1, data = gridobj$grid)$X
        A0 <- cbind(A0,
                    .grid.order = gridobj$grid$.grid.order,
                    .u.order = gridobj$grid$.u.order)
        A1 <- cbind(A1,
                    .grid.order = gridobj$grid$.grid.order,
                    .u.order = gridobj$grid$.u.order)
    } else {
        m0 <- as.formula(paste(gsub("\\s+", " ",
                                    Reduce(paste, deparse(m0))), "+", uname))
        m1 <- as.formula(paste(gsub("\\s+", " ",
                                    Reduce(paste, deparse(m1))), "+", uname))
        ## m0 <- update(m0, as.formula(paste("~ . +", uname)))
        ## m1 <- update(m1, as.formula(paste("~ . +", uname)))
        A0 <- design(formula = m0, data = gridobj$grid)$X
        A1 <- design(formula = m1, data = gridobj$grid)$X
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
        ## Experimenting ------------------------------------------
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
        ## End of experiment --------------------------------------
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
                            ## Experimenting ---------------------
                            isFactor  <- FALSE
                            isBoolean <- FALSE
                            iNamePos  <- NULL
                            ## ## Special case for factors
                            ## if (grepl("factor\\([[:alnum:]]*\\)", iName)) {
                            ##     isFactor <- TRUE
                            ##     iNamePos <-
                            ##         which(gsub(
                            ##             "factor\\([[:alnum:]]*\\)[[:alnum:]]*",
                            ##             "",
                            ##             colnames(get(paste0("A", d)))) == "")
                            ## }
                            ## ## Special case for boolean
                            ## if (!isFactor) {
                            ##     for (op in c("==", "!=", ">", ">=", "<", "<=")){
                            ##         if (grepl(op, iName)) {
                            ##             isBoolean <- TRUE
                            ##             iName <- parenthBoolean(iName)
                            ##             print("NEW INAME BOOL")
                            ##             print(iName)
                            ##             ## iNamePos <-
                            ##             ##     which(colnames(get(paste0("A", d)))
                            ##             ##           == paste0(iName, "TRUE"))
                            ##             break
                            ##         }
                            ##     }
                            ## }
                            ## ## Standard case
                            ## if (!isFactor && !isBoolean) {
                            ##     iNamePos <-
                            ##         which(colnames(get(paste0("A", d))) ==
                            ##               iName)
                            ## }
                            print("iname")
                            print(iName)
                            print("namesA")
                            print(namesA)
                            iNameList <- unlist(strsplit(iName, ":"))
                            namesAscore <- rep(0, times = length(namesA))
                            for (q in iNameList) {
                                if (substr(q, nchar(q), nchar(q)) == ")") {
                                    q <- gsub("\\)", "\\\\)",
                                              gsub("\\(", "\\\\(", q))
                                    namesApos <-
                                        as.integer(
                                            grepl(paste0(q, "[[:alnum:]]*"),
                                                  namesA))

                                } else {
                                    namesApos <-
                                        as.integer(
                                            grepl(q, namesA))
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
                                bmatTmp[, 4:ncol(bmatTmp)] <-
                                    sweep(x = bmatTmp[, 4:ncol(bmatTmp)],
                                          MARGIN = 1,
                                          STATS = bmatTmp[, namesA[r]],
                                          FUN = "*")
                                namesB <- paste0(colnames(bmatTmp)[4:ncol(bmatTmp)],
                                                 ":", namesA[r])
                                colnames(bmatTmp)[4:ncol(bmatTmp)] <- namesB
                                newA <- merge(get(paste0("A", d)),
                                              bmatTmp[, c(".grid.order",
                                                          namesB)],
                                              by = ".grid.order")
                                ## newA <- newA[, c(namesA, namesB)]
                                assign(paste0("A", d), newA)
                            }
                            ## End experimenting -----------------

                            ## Original -------------------------------------
                            ## bmat <-
                            ##     merge(
                            ##         get(paste0("A", d))[, c(uname,
                            ##                                 iName,
                            ##                                 ".grid.order")],
                            ##         bmat, by = uname)
                            ## bmat[, 4:ncol(bmat)] <-
                            ##     sweep(x = bmat[, 4:ncol(bmat)],
                            ##           MARGIN = 1,
                            ##           STATS = bmat[, iName],
                            ##           FUN = "*")
                            ## namesB <- paste0(colnames(bmat)[4:ncol(bmat)],
                            ##                  ":", iName)
                            ## colnames(bmat)[4:ncol(bmat)] <- namesB
                            ## newA <- merge(get(paste0("A", d)),
                            ##               bmat[, c(".grid.order", namesB)],
                            ##               by = ".grid.order")
                            ## newA <- newA[, c(namesA, namesB)]
                            ## assign(paste0("A", d), newA)
                            ## End original -------------------------------------

                        } else {
                            namesA <- colnames(get(paste0("A", d)))
                            namesB <- paste0(colnames(bmat)[2:ncol(bmat)],
                                             ":", iName)
                            colnames(bmat)[2:ncol(bmat)] <- namesB
                            newA <- merge(get(paste0("A", d)),
                                          bmat,
                                          by = uname)
                            ## newA <- newA[, c(namesA, namesB)]
                            assign(paste0("A", d), newA)
                        }
                    }
                }
            }
        }
        rownames(A0) <- A0[, ".grid.order"]
        rownames(A1) <- A1[, ".grid.order"]
    }
    A0 <- A0[order(A0[, ".grid.order"]), ]
    A1 <- A1[order(A1[, ".grid.order"]), ]

    ## Rename columns so they match with the names in vectors gstar0
    ## and gstar1 (the problem stems from the unpredictable ordering
    ## of variables in interaction terms).
    for (d in c(0, 1)) {
        Amat <- get(paste0("A", d))
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
    }

    ## Some columns maybe missing relative to gstar0/gstar1 becuase
    ## the grid is not large enough, and so does not contain all
    ## factor variables. Fill these columns with 0

    missingA0 <- !(names(gstar0) %in% colnames(A0))
    for (i in names(gstar0)[missingA0]) {
        A0[, i] <- 0
    }   
    missingA1 <- !(names(gstar1) %in% colnames(A1))
    for (i in names(gstar1)[missingA1]) {
        A1[, i] <- 0
    }

    print(names(gstar0))
    stop('end of test')
    
    ## keep only the columns that are in the MTRs (A0 and A1 matrices
    ## potentially include extraneous columns)
    A0 <- as.matrix(A0[, names(gstar0)])
    A1 <- as.matrix(A1[, names(gstar1)])
    colnames(A0) <- names(gstar0)
    colnames(A1) <- names(gstar1)

    ## generate null objects
    bdA     <- NULL
    monoA   <- NULL
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
                        "mte.lb", "mte.ub")
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
        monolist  <- c("m0.dec", "m0.inc",
                       "m1.dec", "m1.inc",
                       "mte.dec", "mte.inc")
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

    ## Update bound sequence counts
    if (!is.null(bdA$lb0seq)) lb0seq <- bdA$lb0seq
    if (!is.null(bdA$lb1seq)) lb1seq <- bdA$lb1seq
    if (!is.null(bdA$lbteseq)) lbteseq <- bdA$lbteseq
    if (!is.null(bdA$ub0seq)) ub0seq <- bdA$ub0seq
    if (!is.null(bdA$ub1seq)) ub1seq <- bdA$ub1seq
    if (!is.null(bdA$ubteseq)) ubteseq <- bdA$ubteseq

    output <- combinemonobound(bdA, monoA)
    output$gridobj <- gridobj
    output$lb0seq  <- lb0seq
    output$lb1seq  <- lb1seq
    output$lbteseq <- lbteseq
    output$ub0seq  <- ub0seq
    output$ub1seq  <- ub1seq
    output$ubteseq <- ubteseq

    boundLength <- length(lb0seq) + length(lb1seq) + length(lbteseq) +
        length(ub0seq) + length(ub1seq) + length(ubteseq)

    ## Update monotonicity sequence counts
    if (!is.null(monoA$mono0seq)) {
        mono0seq <- monoA$mono0seq
        mono0seq[, 1] <- mono0seq[, 1] + boundLength
        output$mono0seq <- mono0seq
    }
    if (!is.null(monoA$mono1seq)) {
        mono1seq <- monoA$mono1seq
        mono1seq[, 1] <- mono1seq[, 1] + boundLength
        output$mono1seq <- mono1seq
    }
    if (!is.null(monoA$monoteseq)) {
        monoteseq <- monoA$monoteseq
        monoteseq[, 1] <- monoteseq[, 1] + boundLength
        output$monoteseq <- monoteseq
    }
    return(output)
}
