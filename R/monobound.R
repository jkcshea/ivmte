#' Obtain Tukey half-space quantiles
#'
#' This function calculates the Tukey half-space quantiles for
#' multivariate random variables.
#'
#' @param x vector of values.
#' @param data a matrix of data, characterizing the empirical
#'     distribution.
#' @return scalar, representing the quantile of \code{x} under the
#'     empirical distribution characterized by \code{data}.
tukeydist <- function(x, data) {
    if (is.null(dim(data))) {
        data <- as.matrix(data)
        ineqvec <- mapply("<=", split(data, seq(1, nrow(data))), x)
        F <- sum(ineqvec)/length(ineqvec)
    } else {
        ineqvec <- lapply(split(data, seq(1, nrow(data))), "<=", x)
        ineqvec <- lapply(ineqvec, min)
        F <- Reduce("+", ineqvec) / length(ineqvec)
    }
    return(F)
}

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

    map   <- NULL
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

    ## Generate place holders for the matrices representing monotonicity
    monoA0  <- NULL
    monoA1  <- NULL
    monoAte <- NULL
    mono0z  <- NULL
    mono1z  <- NULL
    monotez <- NULL
    mono0s  <- NULL
    mono1s  <- NULL
    monotes <- NULL
    monoA0seq <- NULL
    monoA1seq <- NULL
    monoAteseq <- NULL
    countseq  <- 0

    monomap <- NULL
    umap <- NULL

    ## This matrix should include all the additions 0s on the left
    ## columns
    sn <- length(sset)
    namesA0 <- colnames(A0)
    namesA1 <- colnames(A1)
    namesA  <- c(seq(1, 2 * sn),
                 namesA0,
                 namesA1)
    ## Generate the constraint matrices ("A" in Gurobi)
    ## corresponding to the monotonicity constraints
    ## A matrix for monotonicity of m0
    if ((hasArg(m0.inc) && m0.inc == TRUE) |
        (hasArg(m0.dec) && m0.dec == TRUE)) {
        monoA0 <- A0[uMaxIndex, ] - A0[uMinIndex, ]
        monoA0 <- cbind(matrix(0, nrow = nrow(monoA0), ncol = 2 * sn),
                        monoA0,
                        matrix(0, nrow = nrow(monoA0), ncol = ncol(A1)))
        mono0z <- replicate(nrow(monoA0), 0)
        monoA0seq <- seq(1, nrow(monoA0))
        colnames(monoA0) <- namesA
        countseq <- countseq + nrow(monoA0)
        monomap <- c(monomap, gridobj$map[uMinIndex])
        umap <- rbind(umap, cbind(gridobj$grid[uMinIndex, uname],
                                  gridobj$grid[uMaxIndex, uname]))
    }
    ## A matrix for monotonicity of m1
    if ((hasArg(m1.inc) && m1.inc == TRUE) |
        (hasArg(m1.dec) && m1.dec == TRUE)) {
        monoA1 <- A1[uMaxIndex, ] - A1[uMinIndex, ]
        monoA1 <- cbind(matrix(0, nrow = nrow(monoA1), ncol = 2 * sn),
                        matrix(0, nrow = nrow(monoA1), ncol = ncol(A1)),
                        monoA1)
        mono1z <- replicate(nrow(monoA1), 0)
        monoA1seq <- seq(1, nrow(monoA1)) + countseq
        colnames(monoA1) <- namesA
        countseq <- countseq + nrow(monoA1)
        monomap <- c(monomap, gridobj$map[uMinIndex])
        umap <- rbind(umap, cbind(gridobj$grid[uMinIndex, uname],
                                  gridobj$grid[uMaxIndex, uname]))
    }
    ## A matrix for monotonicity of m1 - m0
    if ((hasArg(mte.inc) && mte.inc == TRUE) |
        (hasArg(mte.dec) && mte.dec == TRUE)) {
        monoAte0 <- -A0[uMaxIndex, ] - A0[uMinIndex, ]
        monoAte1 <- A1[uMaxIndex, ] - A1[uMinIndex, ]
        monoAte <- cbind(matrix(0, nrow = nrow(monoA0), ncol = 2 * sn),,
                         monoAte0,
                         monoAte1)
        monotez <- replicate(nrow(monoAte), 0)
        colnames(monoAte) <- namesA
        monoAteseq <- seq(1, nrow(monoAte)) + countseq
        monomap <- c(monomap, gridobj$map[uMinIndex])
        umap <- rbind(umap, cbind(gridobj$grid[uMinIndex, uname],
                                  gridobj$grid[uMaxIndex, uname]))
    }
   
    ## Now generate the model sense vectors
    if (try(m0.inc, silent = TRUE) == TRUE) {
        mono0s <- replicate(nrow(monoA0), ">=")
    }
    if (try(m0.dec, silent = TRUE) == TRUE) {
        mono0s <- replicate(nrow(monoA0), "<=")
    }
    if (try(m1.inc, silent = TRUE) == TRUE) {
        mono1s <- replicate(nrow(monoA1), ">=")
    }
    if (try(m1.dec, silent = TRUE) == TRUE) {
        mono1s <- replicate(nrow(monoA1), "<=")
    }
    if (try(mte.inc, silent = TRUE) == TRUE) {
        monotes <- replicate(nrow(monoAte), ">=")
    }
    if (try(mte.dec, silent = TRUE) == TRUE) {
        monotes <- replicate(nrow(monoAte), "<=")
    }
    ## Combine matrices and return
    monoA <- rbind(monoA0, monoA1, monoAte)
    monos   <- c(mono0s, mono1s, monotes)
    monorhs <- c(mono0z, mono1z, monotez)
    
    return(list(A = monoA,
                sense = monos,
                rhs = monorhs,
                map = monomap,
                umap = umap,
                mono0seq = monoA0seq,
                mono1seq = monoA1seq,
                monoteseq = monoAteseq))
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
        m0 <- update(m0, as.formula(paste("~ . +", uname)))
        m1 <- update(m1, as.formula(paste("~ . +", uname)))
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
        for (d in 0:1) {
            if (!is.null(basisList[[d + 1]])) {
                for (j in 1:length(splines[[d + 1]])) {
                    for (v in 1:length(splines[[d + 1]][[j]])) {
                        bmat <- cbind(uvec, basisList[[d + 1]][[j]])
                        colnames(bmat)[1] <- uname
                        iName <- splines[[d + 1]][[j]][v]
                        if (iName != "1") {
                            namesA <- colnames(get(paste0("A", d)))
                            bmat <-
                                merge(
                                    get(paste0("A", d))[, c(uname,
                                                            iName,
                                                            ".grid.order")],
                                    bmat, by = uname)

                            bmat[, 4:ncol(bmat)] <-
                                sweep(x = bmat[, 4:ncol(bmat)],
                                      MARGIN = 1,
                                      STATS = bmat[, iName],
                                      FUN = "*")
                            namesB <- paste0(colnames(bmat)[4:ncol(bmat)],
                                             ":", iName)
                            colnames(bmat)[4:ncol(bmat)] <- namesB
                            newA <- merge(get(paste0("A", d)),
                                          bmat[, c(".grid.order", namesB)],
                                          by = ".grid.order")
                            newA <- newA[, c(namesA, namesB)]
                            assign(paste0("A", d), newA)
                        } else {
                            namesA <- colnames(get(paste0("A", d)))
                            namesB <- paste0(colnames(bmat)[2:ncol(bmat)],
                                             ":", iName)
                            colnames(bmat)[2:ncol(bmat)] <- namesB
                            newA <- merge(get(paste0("A", d)),
                                          bmat,
                                          by = uname)
                            newA <- newA[, c(namesA, namesB)]
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

    ## Update monotonicity sequence counts
    if (!is.null(monoA$mono0seq)) mono0seq <- monoA$mono0seq
    if (!is.null(monoA$mono1seq)) mono1seq <- monoA$mono1seq
    if (!is.null(monoA$monoteseq)) monomteseq <- monoA$monoteseq

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

    output$mono0seq <- mono0seq + boundLength
    output$mono1seq <- mono1seq + boundLength
    output$monomteseq <- monomteseq + boundLength
    return(output)
}
