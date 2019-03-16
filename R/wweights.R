#' Auxiliary function that converts a vector of strings into an
#' expression containing variable names.
#' @param vector Vector of variable names (strings).
#' @return An expression for the list of variable names that are not
#'     strings.
#'
#' @examples
#' ivmte:::unstring(c("a", "b"))
unstring <- function(vector) {
    vector <- parse(text = paste0("c(", paste(vector, collapse = ", "), ")"))
    return(vector)
}

#' Auxiliary function that converts an expression of variable names
#' into a vector of strings.
#' @param vector An expression of a list of variable names.
#' @param substitute Boolean option of whether or not we wish to use
#'     the \code{substitute} command when implementing this
#'     function. Note that this substitutes the argument of the
#'     function. If \code{substitute = FALSE}, then the function will
#'     instead treat the arguments as variables, and substitute in
#'     their values.
#' @param command character, the name of the function defining the
#'     vector or list, e.g. "c", "list", "l". This let's the function
#'     determine how many characters in front to remove.
#' @return A vector of variable names (strings).
#'
#' @examples
#' a <- 4
#' b <- 5
#' ivmte:::restring(c(a, b), substitute = TRUE)
#' ivmte:::restring(c(a, b), substitute = FALSE)
restring <- function(vector, substitute = TRUE, command = "c") {
    startPoint <- nchar(command) + 2
    
    if (substitute == TRUE)  vector <- deparse(substitute(vector))
    if (substitute == FALSE) vector <- deparse(vector)
    vector <- substr(vector, startPoint, nchar(vector) - 1)
    vector <- strsplit(vector, ", ")[[1]]
    return(vector)
}

#' Target weight for ATE
#'
#' Function generates the target weight for the ATE.
#' @param data \code{data.frame} on which the estimation is performed.
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
wate1 <- function(data) {
    return(list(lb = replicate(nrow(data), 0),
                ub = replicate(nrow(data), 1),
                mp = 1))
}

#' Target weight for ATT
#'
#' Function generates the target weight for the ATT.
#' @param data \code{data.frame} on which the estimation is performed.
#' @param expd1 Scalar, the probability that treatment is received.
#' @param propensity Vector of propensity to take up treatment.
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
watt1 <- function(data, expd1, propensity) {
    ## it is assumed that propensity will be a vector
    return(list(lb = replicate(nrow(data), 0),
                ub = propensity,
                mp = 1 / expd1))
}

#' Target weight for ATU
#'
#' Function generates the target weight for the ATT.
#' @param data \code{data.frame} on which the estimation is performed.
#' @param expd0 Scalar, the probability that treatment is not
#'     recieved.
#' @param propensity Vector of propensity to take up treatment.
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
watu1 <- function(data, expd0, propensity) {
    return(list(lb = propensity,
                ub = replicate(nrow(data), 1),
                mp = 1 / expd0))
}

#' Target weight for LATE
#'
#' Function generates the target weight for the LATE, conditioned on a
#' specific value of the covariates.
#' @param data \code{data.frame} on which the estimation is performed.
#' @param from Vector of baseline values for the instruments.
#' @param to Vector of comparison values for the instruments.
#' @param Z An expression for the vector of names of instruments.
#' @param model A \code{lm} or \code{glm} object, or a
#'     \code{data.frame}, which can be used to estimate the propensity
#'     to take up treatment for the specified values of the
#'     instruments.
#' @param X Expression of variable names for the non-excluded
#'     variables the user wishes to condition the LATE on.
#' @param eval.X Vector of values the user wishes to condition the
#'     \code{X} variables on.
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
wlate1 <- function(data, from, to, Z, model, X, eval.X) {

    ## Determine the type of model we are working with (data.frame
    ## vs. glm)
    modclass <- class(model)[1]

    zIsVec <- substr(deparse(Z), 1, 2) == "c("
    if (zIsVec) {
        strinst <- restring(Z, substitute = FALSE)
    } else {
        strinst <- deparse(Z)
    }
    strcovar <- NULL

    if (!is.null(X)) {
        strcovar <- restring(X, substitute = FALSE)
        data[, strcovar] <- t(replicate(nrow(data), eval.X))
    }

    ## Predict propensity scores for 'from' case
    if (length(strinst) == 1) {
        data[, strinst] <- replicate(nrow(data), from)
    } else {
        data[, strinst] <- t(replicate(nrow(data), from))
    }

    if (modclass ==  "lm") {
        bfrom <- predict.lm(model, data)
    }

    if (modclass == "glm") {
        bfrom <- predict.glm(model, data,
                             type = "response")
    }

    ## Predict propensity scores for 'from' case
    if (length(strinst) == 1) {
        data[, strinst] <- replicate(nrow(data), to)
    } else {
        data[, strinst] <- t(replicate(nrow(data), to))
    }

    if (modclass ==  "lm") {
        bto   <- predict.lm(model, data)
    }

    if (modclass == "glm") {
        bto   <- predict.glm(model, data,
                           type = "response")
    }

    ## Predict propensity scores using data.frame model
    if (modclass == "data.frame") {

        cond_from <- mapply(function(a, b) paste(a, "==", b), strinst, from)
        cond_from <- paste(cond_from, collapse = " & ")

        cond_to <- mapply(function(a, b) paste(a, "==", b), strinst, to)
        cond_to <- paste(cond_to, collapse = " & ")

        if (!is.null(X)) {
            condX <- mapply(function(a, b) paste(a, "==", b), strcovar, eval.X)
            condX <- paste(condX, collapse = " & ")
            cond_from <- paste(c(cond_from, condX), collapse = " & ")
            cond_to   <- paste(c(cond_to, condX), collapse = " & ")
        }

        pname <- colnames(model)[(!colnames(model) %in% c(strinst, strcovar))]
        bfrom <- subset(model, eval(parse(text = cond_from)))[, pname]
        bto   <- subset(model, eval(parse(text = cond_to)))[, pname]
    }

    lb <- min(bfrom, bto)
    ub <- max(bfrom, bto)

    ## Ensure the bounds are within 0 and 1
    if (lb < 0) {
        lb <- 0
        warning("Propensity scores below 0 set to 0.", immediate. = TRUE)
    }
    if (ub > 1) {
        ub <- 1
        warning("Propensity scores above 1 set to 1.", immediate. = TRUE)
    }

    return(list(lb = replicate(nrow(data), lb),
                ub = replicate(nrow(data), ub),
                mp =  1 / (ub - lb)))
}

#' Target weight for generalized LATE
#'
#' Function generates the target weight for the generalized LATE,
#' where the user can specify the interval of propensity scores
#' defining the compliers.
#' @param data \code{data.frame} on which the estimation is performed.
#' @param ulb Numeric, lower bound of interval.
#' @param uub Numeric, upper bound of interval.
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
wgenlate1 <- function(data, ulb, uub) {
    return(list(lb = replicate(nrow(data), ulb),
                ub = replicate(nrow(data), uub),
                mp = 1 / (uub - ulb)))
}

#' Generating list of target weight functions
#'
#' This function takes in the user-defined target weight functions and
#' the data set, and generates the weight functions for each
#' observation.
#' @param fun custom weight function defined by the user. Arguments of
#'     the weight function must only be names of variables entering
#'     into the function, and can include the unobserved variable.
#' @param fun.name string, name of function.
#' @param uname the name assigned to the unobserved variable entering
#'     into the MTR.
#' @param data a named vector containing the values of the variables
#'     defining the 'fun', excluding the value of the unobservable
#'     (generated from applying split() to a data.frame).
#' @return The weight function 'fun', where all arguments other than
#'     that of the unobserved variable are fixed according to the
#'     vector 'data'.
genWeight <- function(fun, fun.name, uname, data) {
    wArgList <- formalArgs(fun)
    wArgListOth <- wArgList[wArgList != uname]
    wArgListInput <- paste(paste(wArgListOth, "=", data[, wArgListOth]),
                           collapse = ", ")

    new_call <- paste0(fun.name, "(", uname, " = u, ", wArgListInput, ")")

    outFunction <- function(u) {
        eval(parse(text = new_call))
    }
    return(outFunction)
}


#' Check if custom weights are negations of each other
#'
#' This function checks whether the user-declared weights for treated
#' and control groups are in fact negations of each other. This is
#' problematic for the GMM procedure when accounting for estimation
#' error of the target weights.
#' @param data data set used for estimation. The comparisons are made
#'     only on values in the support of the data set.
#' @param target.knots0 user-defined set of functions defining the
#'     knots associated with splines weights for the control
#'     group. The arguments of the function should consist only of
#'     variable names in \code{data}. If the knot is constant across
#'     all observations, then the user can instead submit the value of
#'     the weight instead of a function.
#' @param target.knots1 user-defined set of functions defining the
#'     knots associated with splines weights for the treated
#'     group. The arguments of the function should be variable names
#'     in \code{data}. If the knot is constant across all
#'     observations, then the user can instead submit the value of the
#'     weight instead of a function.
#' @param target.weight0 user-defined weight function for the control
#'     group defining the target parameter. A list of functions can be
#'     submitted if the weighting function is in fact a spline. The
#'     arguments of the function should be variable names in
#'     \code{data}. If the weight is constant across all observations,
#'     then the user can instead submit the value of the weight
#'     instead of a function.
#' @param target.weight1 user-defined weight function for the treated
#'     group defining the target parameter. A list of functions can be
#'     submitted if the weighting function is in fact a spline. The
#'     arguments of the function should be variable names in
#'     \code{data}. If the weight is constant across all observations,
#'     then the user can instead submit the value of the weight
#'     instead of a function.
#' @param N integer, default set to 20. This is the maxmimum number of
#'     points between treated and control groups to compare and
#'     determine whether or not the weights are indeed negations of
#'     one another. If the data set contains fewer than \code{N}
#'     unique values for a given set of variables, then all those
#'     unique values are used for the comparison.
#' @return boolean. If the weights are negations of each other,
#'     \code{TRUE} is returned.
negationCheck <- function(data, target.knots0, target.knots1,
                          target.weight0, target.weight1, N = 20) {

    negation <- TRUE
    ## Check number of knots
    if (length(target.knots0) == length(target.knots1)) {
        i <- 1
        while (i <= length(target.weight0)) {
            if (i < length(target.weight0)) {
                ## Check if knot arguments are the same
                if (!setequal(formalArgs(target.knots0[[i]]),
                              formalArgs(target.knots1[[i]]))) {
                    negation <- FALSE
                    break
                }

                ## If knots arguments are the same, check
                ## if knot values are the same
                if (!setequal(formalArgs(target.knots0[[i]]), "...")) {
                    wKnotVars <- formalArgs(target.knots0[[i]])
                    tdata <- unique(data[, wKnotVars])
                    if (is.null(dim(tdata))) {
                        tdata <- matrix(tdata, ncol = length(wValVars))
                    }
                    kNCheck <- min(N, nrow(tdata))
                    kNSample <- sample(x = seq(1, nrow(tdata)),
                                       size = kNCheck,
                                       replace = FALSE)

                    kNEval0 <- unlist(
                        lapply(X = split(tdata[kNSample, ],
                                         seq(1, kNCheck)),
                               FUN = funEval,
                               fun = target.knots0[[i]],
                               argnames = wKnotVars))

                    kNEval1 <- unlist(
                        lapply(X = split(tdata[kNSample, ],
                                         seq(1, kNCheck)),
                               FUN = funEval,
                               fun = target.knots1[[i]],
                               argnames = wKnotVars))

                    if (! all(kNEval0 == kNEval1)) {
                        negation <- FALSE
                        break
                    }
                } else {
                    if (target.knots0[[i]](0) !=
                                          target.knots1[[i]](0)) {
                        negation <- FALSE
                        break
                    }
                }
            }

            ## Check if weight arguments are the same
            if (!setequal(formalArgs(target.weight0[[i]]),
                          formalArgs(target.weight1[[i]]))) {
                negation <- FALSE
                break
            }

            ## Check if weights are the same
            if (!setequal(formalArgs(target.weight0[[i]]), "...")) {
                wValVars <- formalArgs(target.weight0[[i]])
                tdata <- unique(data[, wValVars])
                if (is.null(dim(tdata))) {
                    tdata <- matrix(tdata, ncol = length(wValVars))
                }

                kNCheck <- min(20, nrow(tdata))
                kNSample <- sample(x = seq(1, nrow(tdata)),
                                   size = kNCheck,
                                   replace = FALSE)

                kNEval0 <- unlist(
                    lapply(X = split(tdata[kNSample, ],
                                     seq(1, kNCheck)),
                           FUN = funEval,
                           fun = target.weight0[[i]],
                           argnames = wValVars))

                kNEval1 <- unlist(
                    lapply(X = split(tdata[kNSample, ],
                                     seq(1, kNCheck)),
                           FUN = funEval,
                           fun = target.weight1[[i]],
                           argnames = wValVars))

                if (! all(kNEval0 == -kNEval1)) {
                    negation <- FALSE
                    break
                }
            } else {
                if (target.weight0[[i]](0) !=
                                       -target.weight1[[i]](0)) {
                    negation <- FALSE
                    break
                }
            }
            i <- i + 1
        }
    } else {
        negation <- FALSE
    }
    return(negation)
}
