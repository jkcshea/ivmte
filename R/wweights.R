#' Auxiliary function that converts a vector of strings into an
#' expression containing variable names.
#' @param vector Vector of variable names (strings).
#' @return An expression for the list of variable names that are not
#'     strings.
#'
#' @examples
#' \code{> unstring(c("a", "b"))}
#' \code{expression(c(a, b))}.
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
#' @return A vector of variable names (strings).
#'
#' @examples
#' \code{> a <- 4}
#' \code{> b <- 4}
#' \code{> restring(c(a, b), substitute = TRUE)}
#' \code{"a" "b"}
#' #' \code{> restring(c(a, b), substitute = TRUE)}
#' \code{"4" "5"}
restring <- function(vector, substitute = TRUE) {
    if (substitute == TRUE)  vector <- deparse(substitute(vector))
    if (substitute == FALSE) vector <- deparse(vector)
    vector <- substr(vector, 3, nchar(vector) - 1)
    vector <- strsplit(vector, ", ")[[1]]
    return(vector)
}

#' Target weight for ATE
#'
#' Function generates the target weight for the ATE.
#' @param data \code{data.frame} on which the estimation is performed.
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
wate1.mst <- function(data) {
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
watt1.mst <- function(data, expd1, propensity) {
    ## it is assumed that propensity will be a vector
    return(list(lb = replicate(nrow(data), 0),
                ub = propensity,
                mp = 1 / expd1))
}

#' Target weight for ATU
#'
#' Function generates the target weight for the ATT.
#' @param data \code{data.frame} on which the estimation is performed.
#' @param expd1 Scalar, the probability that treatment is not
#'     recieved.
#' @param propensity Vector of propensity to take up treatment.
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
watu1.mst <- function(data, expd0, propensity) {
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
wlate1.mst <- function(data, from, to, Z, model, X, eval.X) {

    ## Determine the type of model we are working with (data.frame
    ## vs. glm)
    modclass <- class(model)[1]
    strinst  <- restring(Z, substitute = FALSE)
    strcovar <- NULL

    if (!is.null(X)) {
        strcovar <- restring(X, substitute = FALSE)
        data[, strcovar] <- t(replicate(nrow(data), eval.X))
    }


    ## Predict propensity scores using lm and glm models
    data[, strinst]  <- t(replicate(nrow(data), from))
    if (modclass ==  "lm") {
        bfrom <- predict.lm(model, data)
    }

    if (modclass == "glm") {
        bfrom <- predict.glm(model, data,
                             type = "response")
    }

    data[, strinst]  <- t(replicate(nrow(data), to))
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
#' @param ulb Numeric, upper bound of interval
#' @param ulb Numeric, lower bound of interval
#' @return The bounds of integration over unobservable \code{u}, as
#'     well as the multiplier in the weight.
wgenlate1.mst <- function(data, ulb, uub) {
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
#'     into the function.
#' @param uname the name assigned to the unobserved variable entering
#'     into the MTR.
#' @param data a vector containing the values of the variables
#'     defining the 'fun', excluding the value of the unobservable.
#' @return The weight function 'fun', where all arguments other than
#'     that of the unobserved variable are fixed according to the
#'     vector 'data'.
genWeight <- function(fun, fun.name,  uname, data) {
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


