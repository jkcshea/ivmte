#' Auxiliary function: generate all permutation orderings
#'
#' This function generates every permutation of the first n natural
#' numbers.
#' @param n integer, the first n natural numbers one wishes to
#'     permute.
#' @return a list of all the permutations of the first n natural
#'     numbers.
permuteN <- function(n){
    if (n == 1) {
        return(matrix(1))
    } else {
        sp <- permuteN(n - 1)
        p <- nrow(sp)
        A <- matrix(nrow = n * p, ncol = n)
        for (i in 1:n){
            A[(i - 1) * p + 1:p, ] <- cbind(i, sp + (sp >= i))
        }
        return(A)
    }
}

#' Auxiliary function: generate all permutations of a vector
#'
#' This function generates every permutation of the elements in a
#' vector.
#' @param vector The vector whose elements are to be permuted.
#' @return a list of all the permutations of \code{vector}.
permute <- function(vector) {
    n <- length(vector)
    permList <- permuteN(n)
    lapply(split(permList, seq(1, nrow(permList))),
            function(x) vector[x])
}

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
piv <- function(Y, X, Z, lmcomponents, weights = NULL) {

    ## project regressors x on image of instruments z
    if (ncol(Z) < ncol(X)) {
        stop(gsub("\\s+", " ",
                  paste0("More regressors than instruments in the following
                  IV-like specification: ", formula, ".")))
    }
    fstage <- if (is.null(weights)) lm.fit(Z, X) else lm.wfit(Z, X, weights)
    Xhat   <- as.matrix(fstage$fitted.values)
    colnames(Xhat) <- colnames(X)

    ## main regression
    if (is.null(weights)) {
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
#' @param order integer, default set to \code{NULL}. This is simply an
#'     index of which IV-like specification the estimate corresponds
#'     to.
#' @return Returns a list containing the matrices of IV-like
#'     specifications for \code{D = 0} and \code{D = 1}; and the
#'     estimates of the IV-like estimands.
#'
#' @examples
#' ivEstimate(formula = ey ~ d | z,
#'            data = dtm,
#'            components = l(d),
#'            treat = d,
#'            list = FALSE)
#' @export
ivEstimate <- function(formula, data, subset, components, treat,
                         list = FALSE, order = NULL) {
    formula <- Formula::as.Formula(formula)
    call <- match.call(expand.dots = FALSE)

    ## Select components
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

    ## Address deparse byte limit
    if (length(components > 1)) {
        components <- gsub("\\s+", " ", Reduce(paste, components))
    }

    ## Address factors whose values are listed, e.g. factor(x)-1,
    ## factor(x)-2. General declaration of factors is dealt with
    ## below.
    components <- gsub(") - ", ")", components)

    ## Covert components into a vector of strings
    stringComp <- (substr(components, 1, 2) == "c(" &
        substr(components, nchar(components), nchar(components)) == ")")
    if (stringComp){
        components   <- substr(components, 3, nchar(components) - 1)
        components   <- strsplit(components, ", ")[[1]]
    }

    ## Some interactions  may need  to be  relabled by  reordering the
    ## order of the interaction
    termsR <- attr(terms(formula), "term.labels")
    failTerms <- which(!components %in% termsR)
    for (fail in failTerms) {
        vars <- strsplit(components[fail], ":")[[1]]
        varsPerm <- permute(vars)
        varsPerm <- unlist(lapply(varsPerm,
                                  function(x) paste(x, collapse = ":")))
        correctPos <- which(varsPerm %in% termsR)
        if (length(correctPos) > 0) {
            components[fail] <- varsPerm[correctPos]
        }
    }
    ## Obtain design matrices
    if (list == TRUE) {
        mf <- design(formula, data)
    } else {
        mf <- eval(modcall(call,
                           newcall = design,
                           keepargs = c("formula", "subset"),
                           newargs = list(data = quote(data))))
    }
    instrumented <- !is.null(mf$Z)
    ## Address factors whose values are not explicitly listed (i.e. if
    ## factor(var1) is povided as a component, this will choose every
    ## variable beginning with factor(var1) in the design matrix
    factorPos <- grep("factor(.)", components)
    if (length(factorPos) > 0) {
        factorVars <- components[factorPos]
        factorPos <- sapply(factorVars,
                            function(x) substr(x, nchar(x), nchar(x)) == ")")
        factorVars <- factorVars[factorPos]
        factorVarsGrep <- sapply(factorVars, function(x) {
            x <- gsub("\\(", "\\\\\\(", x)
            x <- gsub("\\)", "\\\\\\)", x)
            x
        })
        xVars <- colnames(mf$X)
        factorVarsPos <- sapply(factorVarsGrep, grep, x = xVars)
        factorVarsFull <- lapply(factorVarsPos, function(x) xVars[x])
        components <- c(components[! components %in% factorVars],
                        unlist(factorVarsFull))
    }
    ## Deal with boolean expressions (the user will have to fllow a
    ## naming convention, e.g. var1==1TRUE.
    for (op in c("==", "!=",
                 ">", ">=", "<", "<=")) {
        boolPos <- grep(op, components)
        if (length(boolPos) > 0) {
            boolVars <- components[boolPos]
            xVars <- colnames(mf$X)
            boolVarsPos <- sapply(components[boolPos], grep, x = xVars)
            boolVarsFull <- lapply(boolVarsPos, function(x) xVars[x])
            components[boolVarsPos] <- unlist(boolVarsFull)
        }
    }

    ## Ensure components are uniquely declared
    components <- unique(components)

    ## Generate the lmcomponents vector
    lmcomponents <- components
    lmcomponents[lmcomponents == "intercept"] <- "(Intercept)"

    ## Obtain s-weights and the beta-hats
    if (!instrumented) {
        ## For the OLS case, we need to obtain additional design
        ## matrices where we fix treatment.
        if (list == TRUE) {
            data[, colnames(data) == treat] <- 0
            mf0 <- design(formula, data)
            data[, colnames(data) == treat] <- 1
            mf1 <- design(formula, data)
        } else {
            data[, colnames(data) == treat] <- 0
            mf0 <- eval(modcall(call,
                                 newcall = design,
                                 keepargs = c("formula", "subset"),
                                 newargs = list(data = quote(data))))
            data[, colnames(data) == treat] <- 1
            mf1 <- eval(modcall(call,
                                 newcall = design,
                                 keepargs = c("formula", "subset"),
                                 newargs = list(data = quote(data))))
        }
        bhat <- piv(mf$Y, mf$X, mf$X, lmcomponents)
        sweight <- olsj(mf$X, mf0$X, mf1$X, components, treat)
    } else {
        if (length(formula)[2] == 2 &
            dim(mf$X)[2] == dim(mf$Z)[2]) {
            bhat <- piv(mf$Y, mf$X, mf$Z, lmcomponents)
            sweight <- ivj(mf$X, mf$Z, components, treat, order)
        } else if (length(formula)[2] == 2 & dim(mf$X)[2] < dim(mf$Z)[2]) {
            bhat <- piv(mf$Y, mf$X, mf$Z, lmcomponents)
            sweight <- tsls(mf$X, mf$Z, components, treat, order)
        } else if (length(formula)[2] == 2 & dim(mf$X)[2] > dim(mf$Z)[2]) {
            stop(gsub("\\s+", " ",
                      paste0("More regressors than instruments in the following
                      IV-like specification: ", deparse(formula), ".")))
        } else {
            stop(gsub("\\s+", " ",
                      paste0("The following IV-like specification is not
                      properly specified: ", deparse(formula), ".")))
        }
    }

    if (!is.null(sweight$errorfactors)) {
        bhat <- bhat[-which(lmcomponents %in% sweight$errorfactors)]
    }

    return(list(sw0 = sweight$s0, sw1 = sweight$s1, betas = bhat))
}
