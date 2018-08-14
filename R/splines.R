#' Spline basis function of order 1
#'
#' This function is the splines basis function of order 1. This
#' function was coded in accordance to Carl de Boor's set of notes on
#' splines, "B(asic)-Spline Basics".
#' @param x vector, the values at which to evaluate the basis
#'     function.
#' @param knots vector, the internal knots.
#' @param i integer, the basis component to be evaluated.
#' @return scalar.
bX <- function(x, knots, i) {
    lb <- knots[i]
    ub <- knots[i + 1]
    ## as.numeric(sapply(x, function(y) (y >= lb & y < ub)))
    as.numeric(x >= lb & x < ub)
}

#' Generating splines weights
#'
#' This function generates the weights required to construct splines
#' of higher order. This function was coded in accordance to Carl de
#' Boor's set of notes on splines, "B(asic)-Spline Basics".
#' @param x vector, the values at which to evaluate the basis
#'     function.
#' @param knots vector, the internal knots.
#' @param i integer, the basis component to be evaluated.
#' @param order integer, the order of the basis. Do not confuse this
#'     with the degree of the splines, i.e. order = degree + 1. 
#' @return scalar. 
weights <- function(x, knots, i, order) {
  
    if (knots[i] != knots[i + order - 1]) {
        (x - knots[i]) / (knots[i + order - 1] - knots[i])
    } else {
        rep(0, length(x))
    }
}

#' Constructing higher order splines
#'
#' This function recursively constructs the higher order splines
#' basis. Note that the function does not take into consideration the
#' order of the final basis function. The dimensions of the inputs
#' dicate this, and are updated in each iteration of the
#' recursion. The recursion ends once the row number of argument
#' \code{bmat} reaches 1. This function was coded in accordance to
#' Carl de Boor's set of notes on splines, "B(asic)-Spline Basics".
#' @param x vector, the values at which to evaluate the basis
#'     function.
#' @param bmat matrix. Each column of \code{bmat} corresponds to an
#'     element of argument \code{x}. Each row corresponds to the
#'     evaluation of basis component \code{i}, \code{i + 1}, .... The
#'     recursive nature of splines requires that we initially evaluate
#'     the basis functions for components \code{i}, ..., \code{i +
#'     degree of spline}. Each iteration of the recursion reduces the
#'     row of \code{bmat} by 1. The recursion terminates once
#'     \code{bmat} has only a single row.
#' @param knots vector, the internal knots.
#' @param i integer, the basis component of interest.
#' @param current.order integer, the current order associated with the
#'     argument \code{bmat}.
#' @return vector, the evaluation of the spline at each value in
#'     vector \code{x}.
splineUpdate <- function(x, bmat, knots, i, current.order) {
    if (is.null(dim(bmat)) & length(x) > 1) {
        return(bmat)
    }

    if (is.null(dim(bmat)) & length(x) == 1) {
        if (length(bmat) == 1) {
            return(bmat)
        } else {
            bmat <- as.matrix(bmat)
        }
    }

    if (dim(bmat)[1] == 1) return(bmat)
   
    ## print("THIS IS THE SEUQUENCE OF iS I WILL USE")
    ## print(seq(i, (i + nrow(bmat) - 1)))
    ## Update bmat      
    wmat <- as.matrix(sapply(x,
                             function(y) {
                                 sapply(X = seq(i, (i + nrow(bmat) - 1)),
                                        FUN = weights,
                                        x = y,
                                        knots = knots,
                                        order = current.order + 1)}))
    
    bmat1 <- as.matrix(wmat * bmat)
    bmat2 <- as.matrix((1 - wmat) * bmat)
    bmat <- bmat1[-nrow(bmat), ] + bmat2[-1, ]
    
    ## Impose recursion
    splineUpdate(x, bmat, knots, i, current.order + 1)
}

#' Evaluating splines basis functions
#'
#' This function evaluates the splines basis functions. Unlike the
#' \code{bSpline} in the \code{splines2} package, this function
#' returns the value of a single spline basis, rather than a vector of
#' values for all the spline basis functions.
#' @param x vector, the values at which to evaluate the basis
#'     function.
#' @param knots vector, the internal knots.
#' @param degree integer, the degree of the splines.
#' @param intercept boolean, default set to \code{TRUE}. This includes
#'     an additional component to the basis splines so that the
#'     splines are a partition of unity (i.e. the sum of all
#'     components equal to 1).
#' @param i integer, the basis component to be evaluated.
#' @param boundary.knots vector, default is \code{c(0, 1)}.
#' @return scalar.
splinesBasis <- function(x, knots, degree, intercept = TRUE, i,
                            boundary.knots = c(0, 1)) {
    if (i > degree + length(knots) + intercept) {
        return(NA)
    }
    else {   
        knots <- sort(c(knots, boundary.knots[2], rep(boundary.knots, degree)))
        if(intercept == TRUE) knots <- c(boundary.knots[1], knots)
        
        bmat <- as.matrix(sapply(x,
                                 function(y) {
                                     sapply(X = seq(i, i + degree),
                                            FUN = bX,
                                            x = y,
                                            knots = knots)}))

        if (ncol(bmat) != length(x) &
            nrow(bmat) == length(x)) bmat <- t(bmat)

        bmat <- splineUpdate(x, bmat, knots, i = i, current.order = 1)
        return(bmat)
    }
}

#' (Alternative) Defining single splines basis functions, with
#' interactions
#'
#' This function returns a numerically integrable function
#' corresponding to a single splines basis function. It was not
#' implemented because it was slower than using the function from the
#' \code{splines2} package.
#' @param splineslist a list of splines commands and names of
#'     variables that interact with the splines. This is generated
#'     using the command \code{\link{removeSplines}}.
#' @param j the index for the spline for which to generate the basis
#'     functions.
#' @param l the index for the basis.
#' @param v a constant that multiplies the spline basis.
#' @return a vectorized function corresponding to a single splines
#'     basis function that can be numerically integrated.
altDefSplinesBasis <- function(splineslist, j, l, v = 1) {
    cmd <- names(splineslist)[j]
    cmd <- gsub("uSplines\\(", "splinesBasis(x = u, i = l, ", cmd)
    
    fun <- function(u) {
        v * eval(parse(text = cmd))
    }
    return(fun)
}
