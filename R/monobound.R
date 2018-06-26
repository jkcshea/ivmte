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
#' the support of the unobservable variable. A cartesian product of
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
gengrid.mst <- function(index, xsupport, usupport, uname) {

    subsupport <- xsupport[index, ]
    if (is.null(dim(subsupport))) {
        subsupport <- data.frame(subsupport)
        colnames(subsupport) <- colnames(xsupport)
    }
    subsupport$.grid.index <- index

    ## generate a record for which rows correspond to which
    ## index---this will be useful for the audit.
    supportrep <- do.call("rbind",
                          replicate(length(usupport),
                                    subsupport,
                                    simplify = FALSE))
    uvecrep <- rep(usupport, each = length(index))
    grid <- cbind(supportrep, uvecrep, seq(1, length(uvecrep)))
    rownames(grid) <- grid$.grid.order
    map <- grid$.grid.index

    grid$.grid.index <- NULL
    colnames(grid) <- c(colnames(xsupport), uname, ".grid.order")

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
genboundA.mst <- function(A0, A1, sset, gridobj, uname,
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
        }
        if (is.numeric(try(m0.lb, silent = TRUE))) {
            lbdA0 <- bdA0
            m0lb  <- replicate(nrow(A0), m0.lb)
            m0lbs <- replicate(nrow(A0), ">=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
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
        }
        if (is.numeric(try(m1.lb, silent = TRUE))) {
            lbdA1 <- bdA1
            m1lb  <- replicate(nrow(A1), m1.lb)
            m1lbs <- replicate(nrow(A1), ">=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
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
        }
        if (is.numeric(try(mte.lb, silent = TRUE))) {
            lbdAte <- bdAte
            telb  <- replicate(nrow(A1), mte.lb)
            telbs <- replicate(nrow(A1), ">=")
            map <- c(map, gridmap)
            umap <- c(umap, grid[, uname])
        }
    }

    ## Combine matrices and return
    bdA <- rbind(lbdA0,  lbdA1,  lbdAte,
                 ubdA0,  ubdA1,  ubdAte)
    bds   <- c(m0lbs, m1lbs, telbs,
               m0ubs, m1ubs, teubs)
    bdrhs <- c(m0lb,  m1lb,  telb,
               m0ub,  m1ub,  teub)


    return(list(A = bdA, sense = bds, rhs = bdrhs, map = map, umap = umap))
}

#' Taking first differences of constraint matrices
#'
#' This function takes in the matrix of values of the MTR evaluated
#' over the grid generated for the audit procedure. The grid is
#' ordered according to the covariates first, and then by the
#' unobservables (this is done in by \code{\link{genmonoA}}). This
#' function takes the first difference of the unobservables within
#' each set of values for the covariates. This is sufficient to
#' generate the monotoncity constraint matrix.
#' @param A a design matrix that evaluates the MTRs over the grid
#'     generated for the audit procedure.
#' @param monogrid the grid generated for the audit procedure, sorted
#'     by the values of the covariates, and monotone increasing in the
#'     unobservable.
#' @param sn the number of binding constraints in the S-set.
#' @param d indicator for treatment group (\code{d=1}) versus control
#'     group (\code{d = 0}).
#' @param ndcols number of terms in the MTR for the other experimental
#'     group. This is used to generate a matrix that is of the correct
#'     dimension.
#' @return a matrx representing the monotonicity restrictions.
diffA <- function(A, monogrid, sn, d, ndcols) {

    A_mono <- A[rownames(monogrid),]
    A_max  <- A_mono[!as.logical(maxminmatch(monogrid,
                                             ".mst.monoc",
                                             ".mst.monog",
                                             min)), ]
    A_min  <- A_mono[!as.logical(maxminmatch(monogrid,
                                             ".mst.monoc",
                                             ".mst.monog",
                                             max)), ]
    mono   <- A_max - A_min
    if (is.null(dim(mono))) mono <- t(as.matrix(mono))

    if (d == 0) monoA <- cbind(matrix(0, nrow = nrow(mono), ncol = 2 * sn),
                               mono,
                               matrix(0, nrow = nrow(mono), ncol = ndcols))
    if (d == 1) monoA <- cbind(matrix(0, nrow = nrow(mono), ncol = 2 * sn),
                               matrix(0, nrow = nrow(mono), ncol = ndcols),
                               mono)
    return(monoA)
}

#' Stacking monotonicity constraint matrices and vectors
#'
#' This function generates the objects in the LP problem associated
#' with the monotonicity constraints declared by the user. This
#' function simply stacks the matrices corresponding to the
#' monotonicity constraints declared by the user. It also stacks the
#' RHS vector associated with the monotonicity constraints, and stacks
#' the vector of inequalities. It is called by the wrapper function
#' \code{genmonoA.mst}.
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
#' @param monogrid a list containing the grid over which the
#'     monotonicity and boundedness conditions are imposed on.
#' @param gstar0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param gstar1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param m0.dec boolean, set to TRUE if MTR for \code{D = 0} should
#'     be monotone decreasing in the unobservable.
#' @param m0.inc boolean, set to TRUE if MTR for \code{D = 0} should
#'     be monotone increasing in the unobservable.
#' @param m1.dec boolean, set to TRUE if MTR for \code{D = 1} should
#'     be monotone decreasing in the unobservable.
#' @param m1.inc boolean, set to TRUE if MTR for \code{D = 1} should
#'     be monotone increasing in the unobservable.
#' @param mte.dec boolean, set to TRUE if MTE should be monotone
#'     decreasing in the unobservable.
#' @param mte.inc boolean, set to TRUE if MTE should be monotone
#'     decreasing in the unobservable.
#' @return a constraint matrix for the LP problem, the associated
#'     vector of inequalities, and the RHS vector in the inequality
#'     constraint. The objects pertain only to the monotonicity
#'     constraints declared by the user.
stackA.mst <- function(A0, A1, sset, monogrid, gstar0, gstar1,
                         m0.dec, m0.inc, m1.dec, m1.inc, mte.dec,
                         mte.inc) {

    sn <- length(sset)

    namesA0 <- colnames(A0)
    namesA1 <- colnames(A1)
    namesA  <- c(seq(1, 2 * sn),
                 namesA0,
                 namesA1)

    ## Generate place holders for the matriecs representing monotonicity
    monoA0  <- NULL
    monoA1  <- NULL
    monoAte <- NULL
    mono0z  <- NULL
    mono1z  <- NULL
    monotez <- NULL
    mono0s  <- NULL
    mono1s  <- NULL
    monotes <- NULL

    ## Generate the constraint matrices ("A" in Gurobi)
    ## corresponding to the monotonicity consraints
    ## A matrix for monotonicity of m0
    if (hasArg(m0.inc) | hasArg(m0.dec)) {
        monoA0 <- diffA(A0, monogrid, sn, 0, length(gstar1))
        mono0z <- replicate(nrow(monoA0), 0)
        colnames(monoA0) <- namesA
    }
    ## A matrix for monotonicity of m1
    if (hasArg(m1.inc) | hasArg(m1.dec)) {
        monoA1 <- diffA(A1, monogrid, sn, 1, length(gstar0))
        mono1z <- replicate(nrow(monoA1), 0)
        colnames(monoA1) <- namesA
    }
    ## A matrix for monotonicity of m1 - m0
    if (hasArg(mte.inc) | hasArg(mte.dec)) {
        monoAte <- diffA(A1, monogrid, sn, 1, length(gstar0)) -
            diffA(A0, monogrid, sn, 0, length(gstar1))
        monotez <- replicate(nrow(monoAte), 0)
        colnames(monoAte) <- namesA
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

    return(list(A = monoA, sense = monos, rhs = monorhs))
}

#' Generate LP components of the monotonicity constraints
#'
#' This function generates the matrix and vectors associated with the
#' monotonicity constraints declared by the user. It takes in a grid
#' of the covariates on which we define the LP contraints, and then
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
#' @param monov name of variable for which the monotonicity conditions
#'     applies to.
#' @return constraint matrix for the LP problem. The matrix pertains
#'     only to the monotonicity conditions on the MTR and MTE declared
#'     by the user.
genmonoA.mst <- function(A0, A1, sset, gridobj, gstar0, gstar1,
                             m0.dec, m0.inc, m1.dec, m1.inc, mte.dec,
                             mte.inc, monov) {

    ## Order columns in the grid so we can take first differences
    ## below (this is one way to construct the constraint matrices
    ## for monotonicity)
    grid <- gridobj$grid
    othercols <- colnames(grid)[(colnames(grid) != monov) &
                                (colnames(grid) != ".grid.order")]
    colorder  <- c(othercols, monov)
    cmdorder <- paste0("order", "(", paste(colorder, collapse = ", "), ")")
    grid$.grid.index <- gridobj$map
    grid <- grid[with(grid, eval(parse(text = cmdorder))), ]

    ## Now group the rows by the combinations of all other variables
    ## other than the variable we are imposing monotonicity on
    if (length(othercols) > 0) {
        monogrid <- grid[alldup(grid[, othercols]), ]
        monogrid <- groupby(monogrid, othercols)
    } else {
        monogrid <- grid
        monogrid$.mst.monog <- 1
        monogrid$.mst.monoc <- seq(1, nrow(monogrid))
    }

    add.audit.x <- length(unique(monogrid$.mst.monog))
    add.audit.i <- nrow(monogrid) / add.audit.x
    uvec <- monogrid[1:add.audit.i, monov]

    umap <- cbind(uvec[-length(uvec)], uvec[-1])
    umap <- do.call("rbind", rep(list(umap), add.audit.x))

    ## obtain the map (simply need to drop one row from every
    ## group, so just drop the row for which count == 1)
    monomap <- monogrid[monogrid$.mst.monoc > 1, ]$.grid.index
    monogrid$.grid.index <- NULL

    ## Now we can construct the matrices for monotonicity
    arglist  <- c("sset",
                  "gstar0", "gstar1",
                  "m0.dec", "m0.inc", "m1.dec",
                  "m1.inc", "mte.dec", "mte.inc")
    monolist <- c("m0.dec", "m0.inc", "m1.dec",
                  "m1.inc", "mte.dec", "mte.inc")

    call <- match.call(expand.dots = FALSE)
    monoAcall <- modcall(call,
                         newcall = stackA.mst,
                         keepargs = arglist,
                         newargs = list(A0 = quote(A0),
                                        A1 = quote(A1),
                                        monogrid = quote(monogrid)))
    monoA <- eval(monoAcall)

    ## expand the map for monotonicity constraints accordingly.
    monoargs <- length(which(match(monolist, names(call), 0) > 0))
    monomap  <- rep(monomap, times = monoargs)
    umap <- do.call("rbind", rep(list(umap), monoargs))

    return(list(A = monoA$A,
                sense = monoA$sense,
                rhs = monoA$rhs,
                map = monomap,
                umap = umap))
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
    }
    if (!is.null(monoA)) {

        mbA    <- rbind(mbA, monoA$A)
        mbs    <- c(mbs, monoA$sense)
        mbrhs  <- c(mbrhs, monoA$rhs)
        mbmap  <- c(mbmap, monoA$map)
        mbumap <- rbind(mbumap, monoA$umap)
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
#'     to generate the grid preceeding the audit.
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
#'     what the user originally input in \code{\link{mst}}, as the
#'     spline components should have been removed. This formula is
#'     simply a linear combination of all covariates that enter into
#'     the original \code{m0} declared by the user in
#'     \code{\link{mst}}.
#' @param m0 one-sided formula for marginal treatment response
#'     function for the treated group. The formula may differ from
#'     what the user originally input in \code{\link{mst}}, as the
#'     spline components should have been removed. This formula is
#'     simply a linear combination of all covariates that enter into
#'     the original \code{m1} declared by the user in
#'     \code{\link{mst}}.
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
        colnames(grid) <- uname
        grid_index <- rownames(grid)
        gridobj <- list(grid = grid,
                        map  = replicate(length(uvec), 1))
    } else {
        gridobj <- gengrid.mst(grid_index,
                               support,
                               uvec,
                               uname)
    }

    if (is.null(splines[[1]]) & is.null(splines[[2]])) {
        A0 <- design.mst(formula = m0, data = gridobj$grid)$X
        A1 <- design.mst(formula = m1, data = gridobj$grid)$X
    } else {
        m0 <- update(m0, as.formula(paste("~ . +", uname)))
        m1 <- update(m1, as.formula(paste("~ . +", uname)))
        
        A0 <- design.mst(formula = m0, data = gridobj$grid)$X
        A1 <- design.mst(formula = m1, data = gridobj$grid)$X

        A0 <- cbind(A0, .grid.order = seq(1, nrow(A0)))
        A1 <- cbind(A1, .grid.order = seq(1, nrow(A1)))
        
        basisList <- list(genBasisSplines.mst(splines = splines[[1]],
                                              x = uvec,
                                              d = 0),
                          genBasisSplines.mst(splines = splines[[2]],
                                              x = uvec,
                                              d = 1))
        
        for (d in 0:1) {
            if (!is.null(basisList[[d + 1]])) {
                for (j in 1:length(splines[[d + 1]])) {
                    for (l in 1:length(splines[[d + 1]][[j]])) {
                        bmat <- cbind(uvec, basisList[[d + 1]][[j]])
                        colnames(bmat)[1] <- uname
                        iName <- splines[[d + 1]][[j]][l]
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

    A0 <- as.matrix(A0[, names(gstar0)])
    A1 <- as.matrix(A1[, names(gstar1)])

    ## generate null objects
    bdA     <- NULL
    monoA   <- NULL

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
                             newcall = genmonoA.mst,
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
    
    return(combinemonobound(bdA, monoA))
}
