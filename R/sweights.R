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
#' @param components Vector of variable names of which user wants the
#'     S-weights for.
#' @param treat Variable name for the treatment indicator.
#' @return A list of two vectors: one is the weight for D = 0, the
#'     other is the weight for D = 1.
olsj.mst <- function(X, components, treat) {

    ## replace intercept name (since user cannot input
    ## parentheses---they don't use strings)
    colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
    cpos <- which(colnames(X) %in% components)
    X0 <- X
    X0[, which(colnames(X) == treat)] <- 0
    X1 <- X
    X1[, which(colnames(X) == treat)] <- 1

    wvec0 <- solve((1 / nrow(X)) * t(X) %*% X) %*% t(X0)
    wvec0 <- extractcols(t(wvec0), cpos)
    colnames(wvec0)  <- components

    wvec1 <- solve((1 / nrow(X)) * t(X) %*% X) %*% t(X1)
    wvec1 <- extractcols(t(wvec1), cpos)
    colnames(wvec1)  <- components

    return(list(s0 = wvec0, s1 = wvec1))
}

#' Wald weights
#'
#' Function generating the S-weighs for Wald estimand. The function
#' operates by using only a single dummy as the instrument. It is up
#' to the user to construct this dummy and subset the data accordingly
#' to get the correct Wald weights.
#' @param D Matrix, one column being 1s, the other being treatment
#'     indicator D.
#' @param Z Matrix, one column being 1s, the other being the
#'     instrument dummy.
#' @return A list of two vectors: one is the weight for D = 0, the
#'     other is the weight for D = 1.
wald.mst <- function(D, Z) {

    D <- D[, colnames(D) != "(Intercept)"]
    Z <- Z[, colnames(Z) != "(Intercept)"]

    wdt  <- cbind(D, Z)
    pz   <- mean(Z)
    ed0  <- mean(wdt[wdt[, 2] == 0, 1])
    ed1  <- mean(wdt[wdt[, 2] == 1, 1])
    wvec <- as.matrix((as.integer(Z == 1) / pz - as.integer(Z == 0) /
                       (1 - pz)) / (ed1 - ed0))
    return(list(s0 = wvec, s1 = wvec))
}

#' IV weights
#'
#' Function generating the S-weights for OLS estimand, with controls.
#' @param X Matrix of covariates, including the treatment indicator.
#' @param Z Matrix of instruments.
#' @param components Vector of variable names of which user wants the
#'     S-weights for.
#' @param treat Variable name for the treatment indicator.
#' @return A list of two vectors: one is the weight for D = 0, the
#'     other is the weight for D = 1.
ivj.mst <- function(X, Z, components, treat) {

    ## replace intercept name (since user cannot input
    ## parentheses---they don't use strings)
    colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
    cpos <- which(colnames(X) %in% components)
    cposcheck <- which(!components %in% colnames(X))
   
    if (length(cposcheck) > 0) {
        errornames  <- paste(components[cposcheck], collapse = ", ")
        matrixnames <- paste(colnames(X), collapse = ", ")
        emessage <- paste0("The following components are not found in the design
                           matrix: ", errornames, ". The variables included in
                           the design matrix are: ", matrixnames, ". Please
                           select the components from the listed variables in
                           the design matrix.")
        emessage <- gsub("\\s+", " ", emessage)
        stop(emessage)
    }

    ## construct weights
    ezx  <- (1 / nrow(X)) * t(Z) %*% X
    wvec <- solve(ezx) %*% t(Z)
    wvec <- extractcols(t(wvec), cpos)
    colnames(wvec)  <- components

    return(list(s0 = wvec, s1 = wvec))
}

#' TSLS weights, with controls
#'
#' Function generating the S-weights for TSLS estimand, with controls.
#' @param X Matrix of covariates, including the treatment indicator.
#' @param Z Matrix of instruments.
#' @param components Vector of variable names of which user wants the
#'     S-weights for.
#' @param treat Variable name for the treatment indicator.
#' @return A list of two vectors: one is the weight for D = 0, the
#'     other is the weight for D = 1.
tsls.mst <- function(X, Z, components, treat) {
    ## replace intercept name (since user cannot input
    ## parentheses---they don't use strings)
    colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
    cpos <- which(colnames(X) %in% components)

    cposcheck <- which(!components %in% colnames(X))
    if (length(cposcheck) > 0) {
        errornames  <- paste(components[cposcheck], collapse = ", ")
        matrixnames <- paste(colnames(X), collapse = ", ")
        emessage <- paste0("The following components are not found in the design
                           matrix: ", errornames, ". The variables included in
                           the design matrix are: ", matrixnames, ". Please
                           select the components from the listed variables in
                           the design matrix.")
        emessage <- gsub("\\s+", " ", emessage)
        stop(emessage)
    }

    ## construct first stage matrix
    exz <- (1 / nrow(X)) * t(X) %*% Z
    ezz <- (1 / nrow(Z)) * t(Z) %*% Z
    pi  <- exz %*% solve(ezz)

    ## construct weights
    wvec <- solve(pi %*% t(exz)) %*% pi %*% t(Z)
    wvec <- extractcols(t(wvec), cpos)
    colnames(wvec) <- components
    return(list(s0 = wvec, s1 = wvec))
}
