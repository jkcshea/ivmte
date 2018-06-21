#' Obtaining IV-like estimands
#'
#' This function performs TSLS to obtain the estimates for the IV-like
#' estimands.
#' 
#' @param Y the vector of outcomes.
#' @param X the matrix of covariates (includes endogenous and
#'     exogenous covariates).
#' @param Z the matrix of instruments (includes exogenous covariates
#'     in the second stage).
#' @param lmcomponents vector of variable names from the second stage
#'     that we want to include in the S-set of IV-like estimands.
#' @param weights vector of weights.
#' @return vector of select coefficient estimates.
#'
#' @examples
#' piv.mst(Y, X, Z, c(d, x1, x2))
piv.mst <- function(Y, X, Z, lmcomponents, weights = NULL) {
    ## FIX: account for weights

    ## project regressors x on image of instruments z
    if (ncol(Z) < ncol(X)) warning("More regressors than instruments")
    fstage <- if(is.null(weights)) lm.fit(Z, X) else lm.wfit(Z, X, weights)
    Xhat   <- as.matrix(fstage$fitted.values)
    colnames(Xhat) <- colnames(X)
  
    ## main regression
    if(is.null(weights)) {
        fit <- lm.fit(Xhat, Y)
    } else {
        fit <- lm.wfit(Xhat, Y, weights)
    }
    return(fit$coefficients[lmcomponents])
}


#' Obtaining IV-like specifications
#'
#' This function estimates the IV-like estimands, as well as generates
#' the IV-like specifications.
#'
#' @param formula formula to be estimated using OLS/IV.
#' @param data \code{data.frame} with which to perform the estimation.
#' @param subset subset condition with which to perform the estimate.
#' @param components vector of variable names whose coefficients we
#'     want to include in the set of IV-like estimands.
#' @param treat name of treatment indicator variable.
#' @param list logical, set to TRUE if this function is being used to
#'     loop over a list of formulas.
#' @return Returns a list containining the matrices of IV-like
#'     specifications for \code{D = 0} and \code{D = 1}; and the
#'     estimates of the IV-like estimands.
#'
#' @examples
#' sweights.mst(y ~ d + x1 | z,
#'              data = data,
#'              subset = x1 < 3,
#'              components = c(d, x1),
#'              treat = d,
#'              list = FALSE)
#'
#' @export 
sweights.mst <- function(formula, data, subset, components = NULL, treat,
                         list = FALSE) {
    
    formula <- Formula::as.Formula(formula)

    call     <- match.call(expand.dots = FALSE)

    ## obtain design matrices
    if (list == TRUE) {
        mf <- design.mst(formula, data)
    } else {
        mf <- eval(modcall(call,
                           newcall = design.mst,
                           keepargs = c("formula", "subset"),
                           newargs = list(data = quote(data))))
    }
    
    instrumented <- !is.null(mf$Z)

    ## select components
    if (list == TRUE) {      
        if (!is.character(components)) {
            components <- deparse(components)
        }
    } else {
        components <- (unlist(lapply(components, deparse)))
        components <- paste0("c(",
                            paste(components, collapse = ", "),
                            ")")

    }
    
    if (components == "NULL") {
        components <- treat
    } else {
        if (substr(components, 1, 2) == "c(" &
            substr(components, nchar(components), nchar(components)) == ")") {
            components   <- substr(components, 3, nchar(components) - 1)
            lmcomponents <- gsub("intercept", "(Intercept)", components)

            components   <- strsplit(components, ", ")[[1]]
            lmcomponents <- strsplit(lmcomponents, ", ")[[1]]
        } else {
            lmcomponents <- components
        }
    }
    
    ## Obtain s-weights and the beta-hats
    if (!instrumented) {
        ## OLS
        bhat <- piv.mst(mf$Y, mf$X, mf$X, lmcomponents)
        sweight <- olsj.mst(mf$X, components, treat)
    } else {
        ## IV
        if (length(formula)[2] == 2 &
            dim(mf$X)[2] == dim(mf$Z)[2]) {
            bhat <- piv.mst(mf$Y, mf$X, mf$Z, lmcomponents)
            sweight <- ivj.mst(mf$X, mf$Z, components, treat)
        }
        
        ## TSLS
        if (length(formula)[2] == 2 & dim(mf$X)[2] < dim(mf$Z)[2]) {
            bhat <- piv.mst(mf$Y, mf$X, mf$Z, lmcomponents)
            sweight <- tsls.mst(mf$X, mf$Z, components, treat)
        }
    }

    return(list(sw0 = sweight$s0, sw1 = sweight$s1, betas = bhat))
}
