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
#' @param order integer, the counter for which IV-like specification
#'     and component the regression is for.
#' @param excluded boolean, to indicate whether or not the regression
#'     involves excluded variables.
#' @return vector of select coefficient estimates.
piv <- function(Y, X, Z, lmcomponents, weights = NULL, order = NULL,
                excluded = TRUE) {
    ## project regressors x on image of instruments z
    if (ncol(Z) < ncol(X)) {
        stop(gsub("\\s+", " ",
                  paste0("More regressors than instruments in the following
                  IV-like specification: ", formula, ".")), call. = FALSE)
    }
    ## Check that the instruments are different from X. If they are
    ## the same, and 'excluded = TRUE', then the user is submitting an
    ## unclear formula.
    if (length(colnames(X)) == length(colnames(Z)) &&
        all(sort(colnames(X)) == sort(colnames(Z))) && excluded == TRUE) {
        stop(gsub("\\s+", " ",
                  paste0("Please ensure that IV-like specifications are not
                          defined such that the first stage covariates are
                          identical to the second stage covariates. In such
                          specifications, simply provide the second stage
                          regression formula.")),
             call. = FALSE)
    }
    fstage <-
        try(if (is.null(weights)) lm.fit(Z, X) else lm.wfit(Z, X, weights),
            silent = TRUE)
    if (class(fstage) == "try-error") {
        errorMessage <- fstage[[1]]
        errorMessage <- gsub("Error in lm.fit\\(Z, X\\) : ", "", errorMessage)
        errorMessage <- gsub("Error in lm.wfit\\(Z, X, weights\\) : ", "",
                             errorMessage)
        stop(gsub("\\s+", " ",
                  paste0("IV-like specification ",
                         order,
                         " generates the following error:\n ",
                         trimws(errorMessage, which = "right"),
                         ". Check the specification and its
                          corresponding subset condition.")), call. = FALSE)
    }
    nullCheck <- apply(X = fstage$coef, MARGIN = 1,
                       FUN = function(x) any(is.na(x)))
    Xhat <- as.matrix(fstage$fitted.values)
    colnames(Xhat) <- colnames(X)
    ## main regression
    if (is.null(weights)) {
        fit <- lm.fit(Xhat, Y)
    } else {
        fit <- lm.wfit(Xhat, Y, weights)
    }
    return(list(coefficients = fit$coefficients[lmcomponents],
                collinearInst = nullCheck))
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
#' dtm <- ivmte:::gendistMosquito()
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
    componentsSplit <- strsplit(components, ":")
    factorPos <- lapply(componentsSplit, function(x) {
        grep(pattern = "^factor\\(.*\\)$", x = x)
    })
    factorSelect <- which(lapply(factorPos, length) > 0)
    if (length(factorSelect > 0)) {
        factorVars <- componentsSplit[factorSelect]
        factorPos <- factorPos[factorSelect]
        xVars <- strsplit(colnames(mf$X), ":")
        compFactors <- NULL
        for (i in 1:length(factorVars)) {
            includeVec <- unlist(lapply(xVars, length)) ==
                length(factorVars[[i]])
            for (j in 1:length(factorVars[[i]])) {
                tmpName <- factorVars[[i]][j]
                tmpName <- gsub("\\(", "\\\\\\(", tmpName)
                tmpName <- gsub("\\)", "\\\\\\)", tmpName)
                tmpIncVec <- unlist(lapply(xVars,
                                           function(x) max(grepl(tmpName, x))))
                includeVec <- includeVec * tmpIncVec
            }
            compFactors <- c(compFactors,
                             colnames(mf$X)[as.logical(includeVec)])
        }
        compFactors <- unique(compFactors)
        components <- c(components[-factorSelect], compFactors)
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
    if (sum(!components[components != "intercept"] %in% colnames(mf$X))) {
        errorVars <- components[!components %in% colnames(mf$X)]
        errorVars <- errorVars[errorVars != "intercept"]
        stop(gsub("\\s+", " ",
                  paste0("The following components are not included in the
                          IV-like specification ", order, ": ",
                         paste(errorVars, sep = ", "), ".")))
    }
    ## Generate the lmcomponents vector
    lmcomponents <- components
    lmcomponents[lmcomponents == "intercept"] <- "(Intercept)"
    ## Obtain s-weights and the beta-hats
    if (!instrumented) {
        ## For the OLS case, we need to obtain additional design
        ## matrices where we fix treatment.
        if (list == TRUE) {
            data[, colnames(data) == treat] <- 0
            mf0 <- design(formula = formula,
                          data = data,
                          treat = treat,
                          orig.names = colnames(mf$X))
            data[, colnames(data) == treat] <- 1
            mf1 <- design(formula = formula,
                          data = data,
                          treat = treat,
                          orig.names = colnames(mf$X))
        } else {
            data[, colnames(data) == treat] <- 0
            mf0 <- eval(modcall(call,
                                 newcall = design,
                                 keepargs = c("formula", "subset"),
                                newargs = list(data = quote(data),
                                               treat = treat,
                                               orig.names = colnames(mf$X))))
            data[, colnames(data) == treat] <- 1
            mf1 <- eval(modcall(call,
                                newcall = design,
                                keepargs = c("formula", "subset"),
                                newargs = list(data = quote(data),
                                               treat = treat,
                                               orig.names = colnames(mf$X))))
        }
        bhat <- piv(mf$Y, mf$X, mf$X, lmcomponents, order = order,
                    excluded = FALSE)$coef
        if (sum(is.na(bhat))) {
            collinearPos <- which(is.na(bhat))
            mfX <- mf$X[, -collinearPos]
            mf0X <- mf0$X[, -collinearPos]
            mf1X <- mf1$X[, -collinearPos]
            collinearCompPos <- which(components %in% names(bhat)[collinearPos])
            if (length(collinearCompPos) > 0) components <-
                                                  components[-collinearCompPos]
            sweight <- olsj(mfX, mf0X, mf1X, components, treat)
            bhat <- bhat[-collinearPos]
        } else {
            sweight <- olsj(mf$X, mf0$X, mf1$X, components, treat)
        }
    } else {
        if (length(formula)[2] == 2 &
            dim(mf$X)[2] <= dim(mf$Z)[2]) {
            ivfit <- piv(mf$Y, mf$X, mf$Z, lmcomponents, order = order,
                         excluded = TRUE)
            bhat <- ivfit$coef
            collinearInst <- ivfit$collinearInst
            if (sum(is.na(bhat)) | sum(collinearInst)) {
                if (sum(is.na(bhat))) {
                    collinearPos <- which(is.na(bhat))
                    mfX <- mf$X[, -collinearPos]
                    collinearCompPos <-
                        which(components %in% names(bhat)[collinearPos])
                    if (length(collinearCompPos) > 0) {
                        components <- components[-collinearCompPos]
                    }
                } else {
                    mfX <- mf$X
                }
                if (sum(collinearInst)) {
                    mfZ <- mf$Z[, !collinearInst]
                } else {
                    mfZ <- mf$Z
                }
                sweight <- tsls(mfX, mfZ, components, treat, order)
                if (sum(is.na(bhat))) bhat <- bhat[-collinearPos]
            } else {
                sweight <- tsls(mf$X, mf$Z, components, treat, order)
            }
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
