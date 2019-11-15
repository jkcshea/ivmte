#' Auxiliary function: generating basis vectors
#'
#' Auxiliary function to generate standard basis vectors.
#' @param pos The position of the non-zero entry/dimension the basis
#'     vector corresponds to
#' @param length Number of dimensions in total/length of vector.
#' @return Vector containing 1 in a single position, and 0 elsewhere.
genej <- function(pos, length) {
    e <- replicate(length, 0)
    e[pos] <- 1
    return(e)
}

#' Auxiliary function: extracting columns by component names
#'
#' Auxiliary function to extract columns from a matrix based on column
#' names.
#' @param M The matrix to extract from.
#' @param components The vector of variable names.
extractcols <- function(M, components) {
    n <- ncol(M)
    emat <- NULL
    for (pos in components) {
        e <- genej(pos, n)
        emat <- cbind(emat, e)
    }
    return(M %*% emat)
}

#' OLS weights
#'
#' Function generating the S-weights for OLS estimand, with controls.
#' @param X Matrix of covariates, including the treatment indicator.
#' @param X0 Matrix of covariates, once fixing treatment to be 1.
#' @param X1 Matrix of covariates, once fixing treatment to be 0.
#' @param components Vector of variable names of which user wants the
#'     S-weights for.
#' @param treat Variable name for the treatment indicator.
#' @param order integer, default set to \code{NULL}. This is simply an
#'     index of which IV-like specification the estimate corresponds
#'     to.
#' @return A list of two vectors: one is the weight for D = 0, the
#'     other is the weight for D = 1.
olsj <- function(X, X0, X1, components, treat, order = NULL) {
    if(!all(colnames(X) == colnames(X0)) |
       !all(colnames(X) == colnames(X1))) {
        stop(gsub("\\s+", " ",
                  "Column names of design matrices to construct OLS
                   weights do not match. This is most likely because
                   the IV-like specification has the treatment variable
                   inside the factor() command. Try rewriting the
                   IV-like specification without placing the treatment
                   variable inside the factor command."))
    }
    ## replace intercept name (since user cannot input
    ## parentheses---they don't use strings)
    colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
    cpos <- unlist(sapply(components, function(x) which(colnames(X) == x)))
    cposcheck <- which(!components %in% colnames(X))
    if (length(cposcheck) > 0) {
        errornames  <- components[cposcheck]
        matrixnames <- paste(colnames(X), collapse = ", ")
        ## Address the case where a factor variable is declared as a
        ## component, but is missing from the design matrix. This will
        ## be assumed to be due to collinearities. The program will
        ## continue on.

        ## errorfactorpos <- grep("factor(.)", errornames)
        ## if (length(errorfactorpos) > 0) {
        ##     errorfactors <- errornames[errorfactorpos]
        ##     components <- components[! components %in% errorfactors]
        ##     errorfactors <- paste(errorfactors, collapse = ", ")
        ##     errornames   <- errornames[-errorfactorpos]
        ##     allfactors <- paste(grep("factor(.)", colnames(X)),
        ##                         collapse = ", ")
        ##     efmessage <- paste0("The following factor components have been
        ##                    dropped: ", paste(errorfactors, collapse = ", "),
        ##                    ".  This may be due to colinearity in the
        ##                    IV specification, or that the factors were not
        ##                    included in the IV specification. The factor
        ##                    variables included in the design matrix are:",
        ##                    allfactors, ".")
        ##     efmessage <- gsub("\\s+", " ", efmessage)
        ##     warning(efmessage)
        ## } else {
        ##     errorfactors <- NULL
        ## }

        ## NOTE: If the code directly above is commented out, then the
        ## code below will stop the function if there are ANY
        ## variables declared as components missing from the design
        ## matrix. This could be either due to errneous input, or
        ## collinearities.

        ## Now address the case where the non-factor
        ## variables declared as components are missing. If any of
        ## them are missing, then it is not assumed to be
        ## collinearity, and the program is stopped.
        if (length(errornames) > 0) {
            errornames <- paste(errornames,
                                  collapse = ", ")
            if (is.null(order)) {
                emessageIV <- "This may be due to collinearity, or that the
                              variable was never included in the IV-like
                              specification."
            } else {
                emessageIV <- paste0("This may be due to collinearity, or that the
                              variable was never included in IV-like
                              specification ", order, ".")
            }
            emessage <-
                paste0("The following components are not found in the design
                       matrix: ", errornames, ". ", emessageIV, " The variables
                       included in the design matrix are: ", matrixnames,
                       ". Please select the components from the listed variables
                       in the design matrix.")
            emessage <- gsub("\\s+", " ", emessage)
            stop(emessage)
        }
    }
    wvec0 <- solve((1 / nrow(X)) * t(X) %*% X) %*% t(X0)
    wvec0 <- extractcols(t(wvec0), cpos)
    colnames(wvec0) <- components
    wvec1 <- solve((1 / nrow(X)) * t(X) %*% X) %*% t(X1)
    wvec1 <- extractcols(t(wvec1), cpos)
    colnames(wvec1) <- components
    return(list(s0 = wvec0, s1 = wvec1))
}

#' TSLS weights, with controls
#'
#' Function generating the S-weights for TSLS estimand, with controls.
#' @param X Matrix of covariates, including the treatment indicator.
#' @param Z Matrix of instruments.
#' @param Z0 Matrix of instruments, fixing treatment to 0.
#' @param Z1 Matrix of instruments, fixing treatment to 1.
#' @param components Vector of variable names of which user wants the
#'     S-weights for.
#' @param treat Variable name for the treatment indicator.
#' @param order integer, default set to \code{NULL}. This is simply an
#'     index of which IV-like specification the estimate corresponds
#'     to.
#' @return A list of two vectors: one is the weight for D = 0, the
#'     other is the weight for D = 1.
tsls <- function(X, Z, Z0, Z1, components, treat, order = NULL) {
    ## replace intercept name (since user cannot input
    ## parentheses---they don't use strings)
    colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
    cpos <- unlist(sapply(components, function(x) which(colnames(X) == x)))
    cposcheck <- which(!components %in% colnames(X))
    if (length(cposcheck) > 0) {
        errornames  <- components[cposcheck]
        matrixnames <- paste(colnames(X), collapse = ", ")
        ## Address the case where a factor variable is declared as a
        ## component, but is missing from the design matrix. This will
        ## be assumed to be due to collinearities. The program will
        ## continue on.

        ## errorfactorpos <- grep("factor(.)", errornames)
        ## if (length(errorfactorpos) > 0) {
        ##     errorfactors <- errornames[errorfactorpos]
        ##     components <- components[! components %in% errorfactors]
        ##     errorfactors <- paste(errorfactors, collapse = ", ")
        ##     errornames   <- errornames[-errorfactorpos]
        ##     allfactors <- paste(grep("factor(.)", colnames(X)),
        ##                         collapse = ", ")
        ##     efmessage <- paste0("The following factor components have been
        ##                    dropped: ", paste(errorfactors, collapse = ", "),
        ##                    ".  This may be due to colinearity in the
        ##                    IV specification, or that the factors were not
        ##                    included in the IV specification. The factor
        ##                    variables included in the design matrix are:",
        ##                    allfactors, ".")
        ##     efmessage <- gsub("\\s+", " ", efmessage)
        ##     warning(efmessage)
        ## } else {
        ##     errorfactors <- NULL
        ## }

        ## NOTE: If the code directly above is commented out, then the
        ## code below will stop the function if there are ANY
        ## variables declared as components missing from the design
        ## matrix. This could be either due to errneous input, or
        ## collinearities.

        ## Now address the case where the non-factor
        ## variables declared as components are missing. If any of
        ## them are missing, then it is not assumed to be
        ## collinearity, and the program is stopped.
        if (length(errornames) > 0) {
            errornames <- paste(errornames,
                                  collapse = ", ")
            if (is.null(order)) {
                emessageIV <- "This may be due to collinearity, or that the
                              variable was never included in the IV-like
                              specification."
            } else {
                emessageIV <- paste0("This may be due to collinearity, or that the
                              variable was never included in IV-like
                              specification ", order, ".")
            }
            emessage <-
                paste0("The following components are not found in the design
                       matrix: ", errornames, ". ", emessageIV, " The variables
                       included in the design matrix are: ", matrixnames,
                       ". Please select the components from the listed variables
                       in the design matrix.")
            emessage <- gsub("\\s+", " ", emessage)
            stop(emessage)
        }
    }
    ## construct first stage matrix
    exz <- (1 / nrow(X)) * t(X) %*% Z
    ezz <- (1 / nrow(Z)) * t(Z) %*% Z
    pi  <- exz %*% solve(ezz)
    ## construct weights
    wvec0 <- solve(pi %*% t(exz)) %*% pi %*% t(Z0)
    wvec0 <- extractcols(t(wvec0), cpos)
    wvec1 <- solve(pi %*% t(exz)) %*% pi %*% t(Z1)
    wvec1 <- extractcols(t(wvec1), cpos)
    colnames(wvec0) <- components
    colnames(wvec1) <- components
    return(list(s0 = wvec0,
                s1 = wvec1))
}
