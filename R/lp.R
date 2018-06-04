#' Constructing LP problem
#'
#' This function takes in the IV estimates from the set of IV
#' regressions decalred by the user, as well as their corresponding
#' moments of the terms in the MTR. These are then used to construct
#' the components that make up the LP problem. Additional constraint
#' matrix is added using \code{mbA} (\code{mb} stands for
#' "monotonicity/boundedness"); extra model sense is added using
#' \code{mbs}; extra RHS values added using \code{mbrhs}).
#' @param sset List of IV-like estimates and the corresponding gamma
#'     terms.
#' @param mbA Matrix used to define the constraints in the LP problem.
#' @param mbs Vector of model sense/inequalities signs used to define
#'     the constraints in the LP problem.
#' @param mbrhs Vector of constants used to define the constraints in
#'     the LP problem.
#' @return A list of matrices and vectors necessary to define an LP
#'     problem for Gurobi.
#'
#' @export
lpsetup.mst <- function(sset, mbA = NULL, mbs = NULL, mbrhs = NULL) {
    
    ## determine lengths
    sn  <- length(sset)
    gn0 <- length(sset$s1$g0)
    gn1 <- length(sset$s1$g1)

    ## generate all vectors/matrices for LP optimizaiton to minimize
    ## observational equivalence
    obj <- c(replicate(sn * 2, 1),
             replicate(gn0 + gn1, 0))
    rhs <- unlist(lapply(sset, function(x) x[["beta"]]))
    sense <- replicate(sn, "=")

    A <- NULL
    scount <- 0
    for (s in names(sset)) {
        avec <- replicate(2 * sn, 0)
        avec[ (2 * scount + 1) :(2 * scount + 2)] <- c(-1, 1) ## (-1 is for w+,
                                                              ##   1 is for w-)
        avec <- c(avec, sset[[s]]$g0, sset[[s]]$g1) ## order of variables determined here
        A <- rbind(A, avec)
        scount <- scount + 1
    }

    ## Add in additional constraints if included
    A     <- rbind(A, mbA)
    sense <- c(sense, mbs)
    rhs   <- c(rhs, mbrhs)

    ## define bounds
    ub <- replicate(ncol(A), Inf)
    lb <- c(replicate(sn * 2, 0), replicate(gn0 + gn1, -Inf))

    return(list(obj = obj,
                rhs = rhs,
                sense = sense,
                A = A,
                ub = ub,
                lb = lb,
                sn = sn,
                gn0 = gn0,
                gn1 = gn1))
}

#' Minimizing violation of observational equivalence
#'
#' Given a set of IV-like estimates and the set of matrices/vectors
#' defining an LP problem, this function minimizes the violation of
#' observational equivalence under the L1 norm.
#' @param sset A list of IV-like estimates and the corresponding gamma
#'     terms.
#' @param lpobj A list of matrices and vectors defining an LP problem.
#' @return A list including the minimum violation of observational
#'     equivalence, the solution to the LP problem, and the status of
#'     the solution.
#'
#' @export
obseqmin.mst <- function(sset, lpobj) {
    
    ## define model
    model <- list()
    model$modelsense <- "min"
    model$obj   <- lpobj$obj 
    model$A     <- lpobj$A
    model$rhs   <- lpobj$rhs
    model$sense <- lpobj$sense
    model$ub    <- lpobj$ub
    model$lb    <- lpobj$lb

    ## solve for min
    result <- gurobi::gurobi(model, list(outputflag = 0))
    min <- result$objval

    ## provide nicer output
    g0sol <- result$x[(2 * lpobj$sn + 1) : (2 * lpobj$sn + lpobj$gn0)]
    g1sol <- result$x[(2 * lpobj$sn + lpobj$gn0 + 1) : (2 * lpobj$sn + lpobj$gn0 + lpobj$gn1)]
    names(g0sol) <- names(sset$gstar$g0)
    names(g1sol) <- names(sset$gstar$g1)
    
    ## return output
    return(list(obj = result$objval,
                g0 = g0sol,
                g1 = g1sol, 
                status = result$status,
                result = result))
}

#' Obtaining TE bounds
#'
#' This function estimates the bounds on the target treatment effect.
#' @param g0 set of expectations for each terms of the MTR for the
#'     control group.
#' @param g1 set of expectations for each terms of the MTR for the
#'     control group.
#' @param sset a list containing the point estimates and gamma
#'     components associated with each element in the S-set.
#' @param lpobj A list of matrices and vectors defining an LP problem.
#' @param threshold tolerance level for how much more the solution is
#'     permitted to violate observational equivalence of the IV-like
#'     estimands.
#' @return a list containing the bounds on the treatment effect; the
#'     coefficients on each term in the MTR associated with the upper
#'     and lower bounds, for both counterfactuals; the optimization
#'     status to the maximization and minimization problems; the LP
#'     problem that the optimizer solved.
#' 
#' @export
bound.mst <- function(g0, g1, sset, lpobj, threshold) {
   
    ## define model
    model <- list()
    model$obj <- c(replicate(2 * lpobj$sn, 0), g0, g1)
    model$rhs <- c(threshold, lpobj$rhs)
    model$sense <- c("<=", lpobj$sense)
    model$ub    <- lpobj$ub
    model$lb    <- lpobj$lb

    avec <- c(replicate(2 * lpobj$sn, 1), replicate(lpobj$gn0 + lpobj$gn1, 0))
    model$A <- rbind(avec, lpobj$A)

    ## find lower bound
    model$modelsense <- "min"
    minresult <- gurobi::gurobi(model, list(outputflag = 0))
    min <- minresult$objval
    ming0 <- minresult$x[(2 * lpobj$sn + 1) : (2 * lpobj$sn + lpobj$gn0)]
    ming1 <- minresult$x[(2 * lpobj$sn + lpobj$gn0 + 1) : (2 * lpobj$sn + lpobj$gn0 + lpobj$gn1)]
    names(ming0) <- names(sset$gstar$g0)
    names(ming1) <- names(sset$gstar$g1)
    
    ## find upper bound
    model$modelsense <- "max"
    maxresult <- gurobi::gurobi(model, list(outputflag = 0))
    max <- maxresult$objval
    maxg0 <- maxresult$x[(2 * lpobj$sn + 1) : (2 * lpobj$sn + lpobj$gn0)]
    maxg1 <- maxresult$x[(2 * lpobj$sn + lpobj$gn0 + 1) : (2 * lpobj$sn + lpobj$gn0 + lpobj$gn1)]
    names(maxg0) <- names(sset$gstar$g0)
    names(maxg1) <- names(sset$gstar$g1)
    
    cat("Min status:", minresult$status, "\n")
    cat("Max status:", maxresult$status, "\n")
    cat("Bound: (", min, ",", max, ")\n")
    
    return(list(max = max,
                maxg0 = maxg0,
                maxg1 = maxg1,
                maxresult = maxresult,
                min = min,
                ming0 = ming0,
                ming1 = ming1,
                minresult = minresult,
                model = model))                
}
