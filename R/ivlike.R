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
#'     that we want to include in the S-set of IV-like estimands. If
#'     \code{NULL} is submitted, then all components will be included.
#' @param weights vector of weights.
#' @param order integer, the counter for which IV-like specification
#'     and component the regression is for.
#' @param excluded boolean, to indicate whether or not the regression
#'     involves excluded variables.
#' @return vector of select coefficient estimates.
piv <- function(Y, X, Z, lmcomponents = NULL, weights = NULL, order = NULL,
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
    if (is.null(dim(fstage$coef))) {
        nullCheck <- is.na(fstage$coef)
    } else {
        nullCheck <- apply(X = fstage$coef, MARGIN = 1,
                           FUN = function(x) any(is.na(x)))
    }
    Xhat <- as.matrix(fstage$fitted.values)
    colnames(Xhat) <- colnames(X)
    ## main regression
    if (is.null(weights)) {
        fit <- lm.fit(Xhat, Y)
    } else {
        fit <- lm.wfit(Xhat, Y, weights)
    }
    if (!is.null(lmcomponents)) {
        return(list(coefficients = fit$coefficients[lmcomponents],
                    collinearInst = nullCheck))
    } else {
        return(list(coefficients = fit$coefficients,
                    collinearInst = nullCheck))
    }
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
            components <- unlist(deparse(components))
        }
    } else {
        components <- unlist(lapply(components, deparse))
        components <- paste0("c(",
                            paste(components, collapse = ", "),
                            ")")
    }
    if (length(components > 1)) {
        components <- gsub("\\s+", " ", Reduce(paste, components))
    }
    ## ------------------------------------------------------
    ## NEW TESTING---COMPLETE REWRITE

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
    xVars <- colnames(mf$X)

    ## Covert components into a vector of strings
    stringComp <- (substr(components, 1, 2) == "c(" &
        substr(components, nchar(components), nchar(components)) == ")")
    if (stringComp){
        components   <- substr(components, 3, nchar(components) - 1)
        components   <- strsplit(components, ", ")[[1]]
    }

    ## Remove dash defining specific factor components
    components <- gsub(") - ", ")", components)
    origComponents <- components
    ## Break apart each component
    componentsSplit <- strsplit(components, ":")

    ## Prepare expansion of components
    componentsFull <- NULL
    componentsType <- NULL
    componentsIndex <- NULL
    ci <- 1
    for (comp in componentsSplit) {
        complist <- list()
        comptype <- NULL
        for (subcomp in comp) {
            ## Exand the undeclared factor variables
            fcheck <- grepl(pattern = "^factor\\(.*\\)$", x = subcomp)
            if (fcheck) {
                fvar <- all.vars(as.formula(paste("~", subcomp)))
                if (length(fvar) > 1) {
                    stop(gsub("\\s+", " ",
                              paste0("Factor expression ",
                                     subcomp, " includes more than
                                     one variable.")))
                }
                grepSubcomp <- gsub("\\(", "\\\\\\(", subcomp)
                grepSubcomp <- gsub("\\)", "\\\\\\)", grepSubcomp)
                grepSubcomp <- gsub("\\.", "\\\\\\.", grepSubcomp)
                fexpanded <- sapply(xVars, function(x) {
                    fpos <- grep(grepSubcomp, strsplit(x, ":"))
                    unlist(strsplit(x, ":"))[fpos]
                })
                fexpanded <- unlist(fexpanded)
                complist[[length(complist) + 1]] <- fexpanded
                comptype <- c(comptype, "f")
            }
            ## Expand the boolean variables
            bcheck <- any(sapply(c("==", "!=", ">", ">=", "<", "<="),
                                 grepl, x = subcomp))
            if (bcheck) {
                bexpanded <- paste0(subcomp, c("TRUE", "FALSE"))
                complist[[length(complist) + 1]] <- bexpanded
                comptype <- c(comptype, "b")
            }
            if (!fcheck && !bcheck) {
                complist[[length(complist) + 1]] <- subcomp
                comptype <- c(comptype, "s")
            }
        }
        ## Create a dictionary:
        ## You need to indicate variables which involve NO expansion at all, SOME expansion, and ALL expansion.
        ## NO expansion variables can never be dropped
        ## SOME expansion variables will be dropped with a warning.
        ## ALL expansions will be dropped without warning.
        comptype <- unique(comptype)
        if (length(comptype) == 1 && comptype == "s") {
            comptype <- 1 ## Specified component---cannot be dropped
        } else if ("s" %in% comptype && ("f" %in% comptype | "b" %in% comptype)) {
            comptype <- 2 ## Contains specified component, with some
                          ## kind of expansion. This will be dropped
                          ## without warning since the specified
                          ## component will appear in every expansion.
        } else if (! "s" %in% comptype) {
            comptype <- 3 ## Only expanded components. Can also be
                          ## dropped without warning.
        }

        ## Construct the full set of components
        if (length(complist) == 1) {
            componentsExp <- complist[[1]]
        } else {
            componentsExp <- complist[[1]]
            for (i in 2:(length(complist))) {
                componentsExp <- outer(componentsExp,
                                       complist[[i]],
                                       paste,
                                       sep = ":")
            }
            componentsExp <- c(componentsExp)
        }
        componentsFull <- c(componentsFull, componentsExp)
        componentsType <- c(componentsType,
                            rep(comptype, length(componentsExp)))
        componentsIndex <- c(componentsIndex,
                             rep(ci, length(componentsExp)))
        ci <- ci + 1
    }
    components <- componentsFull
    ## Some interactions  may need  to be  relabled by  reordering the
    ## order of the interaction
    failTerms <- which(!components %in% xVars)
    for (fail in failTerms) {
        vars <- strsplit(components[fail], ":")[[1]]
        varsPerm <- permute(vars)
        varsPerm <- unlist(lapply(varsPerm,
                                  function(x) paste(x, collapse = ":")))
        correctPos <- which(varsPerm %in% xVars)
        if (length(correctPos) > 0) {
            components[fail] <- varsPerm[correctPos]
        }
    }
    ## Now restrict components
    failTerms <- which(!components %in% xVars)
    if ("intercept" %in% components) {
        failTerms <- failTerms[failTerms != which(components == "intercept")]
    }
    if (length(failTerms) > 0) {
        if (!all(sort(unique(componentsIndex)) ==
                 sort(unique(componentsIndex[-failTerms])))) {
            missingComp <- which(! sort(unique(componentsIndex)) %in%
                                 sort(unique(componentsIndex[-failTerms])))
            missingComp <- origComponents[missingComp]
            stop(gsub("\\s+", " ",
                      paste0("The following components are not included
                          in the IV-like specification: ",
                          paste(missingComp, collapse = ", "),
                          ". Remove them from the components list, or
                          update the corresponding IV-like specification.")))
        }
        components <- components[-failTerms]
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
            collinearVars <- lmcomponents[collinearPos]
            mfX <- mf$X[, -which(colnames(mf$X) %in% collinearVars)]
            mf0X <- mf0$X[, -which(colnames(mf0$X) %in% collinearVars)]
            mf1X <- mf1$X[, -which(colnames(mf0$X) %in% collinearVars)]
            collinearCompPos <- which(components %in% names(bhat)[collinearPos])
            if (length(collinearCompPos) > 0) components <-
                                                  components[-collinearCompPos]
            sweight <- olsj(mfX, mf0X, mf1X, components, treat, order)
            bhat <- bhat[-collinearPos]
        } else {
            sweight <- olsj(mf$X, mf0$X, mf1$X, components, treat, order)
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
                    collinearVars <- lmcomponents[collinearPos]
                    mfX <- mf$X[, -which(colnames(mf$X) %in% collinearVars)]
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
