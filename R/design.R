#' Generating design matrices
#'
#' This function generates the design matrix given an IV
#' specification.
#'
#' @param formula Formula with which to generate the design matrix.
#' @param data \code{data.frame} with which to generate the design
#'     matrix.
#' @param subset Condition to select subset of data.
#' @return Three matrices are returned: one for the outcome variable,
#'     Y; one for the second stage covariates, X; and one for the
#'     first stage covariates, Z.
#'
#' @examples
#' design(formula = ey ~ d | z,
#'            data = dtm,
#'            subset = z %in% c(1, 2))
#'
#' @export
design <- function(formula, data, subset) {
    
    ## Set up model.frame() call
    if (missing(data)) data <- environment(formula)

    ## Convert formula to Formula
    formula <- Formula::as.Formula(formula)
    onesided <- FALSE

    if (length(formula)[1] == 0L) onesided <- TRUE

    call <- match.call()
    design_call <- modcall(call,
                           newcall = model.frame,
                           keepargs = c("data",
                                        "subset"),
                           newargs = list(formula = formula,
                                          drop.unused.levels = TRUE))
    
    mf <- eval(design_call, parent.frame())
                               
    ## extract response, terms, model matrices
    if (onesided == FALSE) Y <- model.response(mf, "numeric")
    if (onesided == TRUE)  Y <- NULL
    mt  <- terms(formula, data = data)
    mtX <- terms(formula, data = data, rhs = 1)
    X   <- model.matrix(mtX, mf)

    if(length(formula)[2] < 2L) {
        mtZ <- NULL
        Z   <- NULL
    } else {
        mtZ <- delete.response(terms(formula, data = data, rhs = 2))
        Z   <- model.matrix(mtZ, mf)
    }

    return(list(Y = Y, X = X, Z = Z))
}
