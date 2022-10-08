#' Generating design matrices
#'
#' This function generates the design matrix given an IV
#' specification.
#'
#' @param formula Formula with which to generate the design matrix.
#' @param data \code{data.frame} with which to generate the design
#'     matrix.
#' @param subset Condition to select subset of data.
#' @param treat The name of the treatment variable. This should only
#'     be passed when constructing OLS weights.
#' @param orig.names character vector of the terms in the final design
#'     matrix. This is required when the user declares an IV-like
#'     formula where the treatment variable is passed into the
#'     \code{factor} function. Since the treatment variable has to be
#'     fixed to 0 or 1, the design matrix will be unable to construct
#'     the contrasts. The argument \code{orig.names} is a vector of
#'     the terms in the IV-like specification prior to fixing the
#'     treatment variable.
#' @return Three matrices are returned: one for the outcome variable,
#'     Y; one for the second stage covariates, X; and one for the
#'     first stage covariates, Z.
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#' design(formula = ey ~ d | z,
#'            data = dtm,
#'            subset = z %in% c(1, 2))
#'
#' @export
design <- function(formula, data, subset, treat, orig.names) {
    ## Set up model.frame() call
    if (missing(data)) data <- environment(formula)
    if (hasArg(treat) && hasArg(orig.names)) {
        ## Formulas where the treatment variable is passed inside of
        ## factor() are problematic when constructing OLS weights (no
        ## contrasts once fixing). The code below converts the
        ## factored treatment variables into booleans, constructs the
        ## matrix, then reverts them back to factors
        origPos <- grepl(paste0("factor\\(", treat, "\\)"), orig.names)
        if (any(origPos)) {
            degenVal <- unique(data[, treat])
            if (length(degenVal) == 1) {
                degenAlt <- 1 - degenVal
                origInter <- orig.names[origPos]
                origInter <- unlist(strsplit(origInter, ";"))
                origInter <- sapply(origInter, function(x) {
                    gsub(paste0("factor\\(", treat, "\\)", degenVal),
                         paste0(treat, " == ", degenVal, "TRUE"), x)
                })
                origInter <- sapply(origInter, function(x) {
                    gsub(paste0("factor\\(", treat, "\\)", degenAlt),
                         paste0(treat, " == ", degenVal, "FALSE"), x)
                })
                formulaStr <- deparse(formula)
                formulaStr <- gsub(paste0("factor\\(", treat, "\\)"),
                                   paste0("(", treat, " == ", degenVal, ")"),
                                   formulaStr)
                formula <- as.formula(formulaStr)
            }
        }
    }
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
    ## Extract response, terms, model matrices
    if (onesided == FALSE) Y <- model.response(mf, "numeric")
    if (onesided == TRUE)  Y <- NULL
    mtX <- terms(formula, data = data, rhs = 1)
    X <- try(model.matrix(mtX, mf), silent = TRUE)
    if (length(X) == 1 && class(X) == "try-error") {
        if (grepl("contrasts", X)) {
            stop(gsub("\\s+", " ",
                      "R was unable to construct a design matrix since the
                      formula involved a contrast for a variable, but the data
                      provided for that variable has no variation.
                      This can occur when constructing the audit grids, which
                      are constructed from a random subset of the data.
                      This error can be resolved by adjusting the seed.
                      A better solution is simply to increase the size of
                      'initgrid.nx' or 'audit.nx'. If custom
                      grids are passed, then make sure the custom grids contain
                      variation for variables declared to be factors in the
                      MTRs."),
                 call. = FALSE)
        } else {
            stop("Unknown error in construction of design matrix.")
        }
    }
    if (hasArg(treat) && hasArg(orig.names)) {
        if (any(origPos)) {
            for (i in 1:length(origInter)) {
                pos <- which(colnames(X) == origInter[i])
                colnames(X)[pos] <- names(origInter)[i]
            }
            if (any(!orig.names %in% colnames(X))) {
                for (i in orig.names[! orig.names %in% colnames(X)]) {
                    tmpXnames <- colnames(X)
                    X <- cbind(X, 0)
                    colnames(X) <- c(tmpXnames, i)
                }
            }
            X <- X[, orig.names]
        }
    }
    if(length(formula)[2] < 2L) {
        mtZ <- NULL
        Z   <- NULL
    } else {
        mtZ <- delete.response(terms(formula, data = data, rhs = 2))
        Z   <- model.matrix(mtZ, mf)
    }
    ## Update column names
    colnames(X) <- parenthBoolean(colnames(X))
    colnames(Z) <- parenthBoolean(colnames(Z))
    return(list(Y = Y, X = X, Z = Z))
}
