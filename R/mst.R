ivmte <- function(bootstraps = 0,
                  ivlike, data, subset, components, propensity, link,
                  treat, m0, m1, uname = u, target, target.weight0,
                  target.weight1, target.knots0 = NULL, target.knots1 = NULL,
                  late.Z, late.from, late.to, late.X,
                  eval.X, genlate.lb, genlate.ub, obseq.tol = 1,
                  grid.nu = 20, grid.nx = 20, audit.nx = 2,
                  audit.nu = 3, audit.max = 5, audit.tol = 1e-08, m1.ub,
                  m0.ub, m1.lb, m0.lb, mte.ub, mte.lb, m0.dec, m0.inc,
                  m1.dec, m1.inc, mte.dec, mte.inc, lpsolver = NULL,
                  point = FALSE, point.itermax = 2,
                  point.tol = 1e-08, noisy = TRUE) {
    ## Match call arguments
    call <- match.call(expand.dots = FALSE)

    estimateCall <- modcall(call,
                            newcall = ivmte.estimate,
                            dropargs = c("bootstraps", "data"),
                            newargs = list(data = quote(data)))

    if (bootstraps == 0) {
        return(eval(estimateCall))
    } else if (point == TRUE & (bootstraps > 0) & (bootstraps %% 1 == 0)) {

        origEstimate <- eval(estimateCall)

        teEstimates  <- NULL
        mtrEstimates <- NULL
        propEstimates <- NULL

        b <- 1
        bootFailN <- 0
        bootFailNote <- ""
        bootFailIndex <- NULL

        while (b <= bootstraps) {
            bdata <- data[sample(seq(1, nrow(data)),
                                 size = nrow(data),
                                 replace = TRUE), ]
            bootCall <- modcall(call,
                                newcall = ivmte.estimate,
                                dropargs = c("bootstraps", "data", "noisy"),
                                newargs = list(data = quote(bdata),
                                               noisy = FALSE))

            bootEstimate <- try(eval(bootCall), silent = TRUE)
            if (is.list(bootEstimate)) {
                teEstimates  <- c(teEstimates, bootEstimate$te)
                mtrEstimates <- cbind(mtrEstimates, bootEstimate$mtr.coef)
                propEstimates <- cbind(propEstimates,
                                       bootEstimate$propensity$model$coef)

                if (noisy == TRUE) {
                    message(paste0("Bootstrap iteration ", b, bootFailNote,
                                   "..."))
                }
                b <- b + 1
                bootFailN <- 0
                bootFailNote <- ""
            } else {
                if (noisy == TRUE) {
                    message(paste0("Bootstrap iteration ", b, bootFailNote,
                                   " error, resampling..."))
                }
                bootFailN <- bootFailN + 1
                bootFailIndex <- unique(c(bootFailIndex, b))
                bootFailNote <- paste0(": resample ",
                                       bootFailN)

            }
        }

        if (length(bootFailIndex) > 0) {
            warning(gsub("\\s+", " ",
                         paste0("Bootstrap iteration(s) ",
                                paste(bootFailIndex, collapse = ", "),
                                " failed. Failed bootstraps are repeated.")))
        }

        bootSE <- sd(teEstimates)
        mtrSE  <- apply(mtrEstimates, 1, sd)
        propSE  <- apply(propEstimates, 1, sd)

        ## Conf. int. 1: quantile method (same as percentile method)
        ci190 <- quantile(x = teEstimates, probs = c(0.05, 0.95), type = 1)
        ci195 <- quantile(x = teEstimates, probs = c(0.025, 0.975), type = 1)

        mtrci190 <- apply(mtrEstimates, 1, quantile, probs = c(0.05, 0.95),
                          type = 1)
        mtrci195 <- apply(mtrEstimates, 1, quantile, probs = c(0.025, 0.975),
                          type = 1)

        propci190 <- apply(propEstimates, 1, quantile, probs = c(0.05, 0.95),
                          type = 1)
        propci195 <- apply(propEstimates, 1, quantile, probs = c(0.025, 0.975),
                          type = 1)

        ## Conf. int. 2: percentile, using Z statistics
        ci290 <- origEstimate$te + c(qnorm(0.05), qnorm(0.95)) * bootSE
        ci295 <- origEstimate$te + c(qnorm(0.025), qnorm(0.975)) * bootSE
        names(ci290) <- c("5%", "95%")
        names(ci295) <- c("2.5%", "97.5%")

        mtrci290 <- sweep(x = tcrossprod(c(qnorm(0.05), qnorm(0.95)), mtrSE),
                          MARGIN = 2, origEstimate$te, FUN = "+")
        mtrci295 <- sweep(x = tcrossprod(c(qnorm(0.025), qnorm(0.975)), mtrSE),
                          MARGIN = 2, origEstimate$te, FUN = "+")

        propci290 <- sweep(x = tcrossprod(c(qnorm(0.05), qnorm(0.95)), propSE),
                          MARGIN = 2, origEstimate$prop$model$coef, FUN = "+")
        propci295 <- sweep(x = tcrossprod(c(qnorm(0.025), qnorm(0.975)), propSE),
                          MARGIN = 2, origEstimate$prop$model$coef, FUN = "+")

        colnames(mtrci290) <- colnames(mtrci190)
        rownames(mtrci290) <- rownames(mtrci190)

        colnames(mtrci295) <- colnames(mtrci195)
        rownames(mtrci295) <- rownames(mtrci195)

        colnames(propci290) <- colnames(propci190)
        rownames(propci290) <- rownames(propci190)

        colnames(propci295) <- colnames(propci195)
        rownames(propci295) <- rownames(propci195)
      
        return(c(origEstimate,
                 list(te.se  = bootSE,
                      mtr.se = mtrSE,
                      prop.se = propSE,
                      te.bootstraps = teEstimates,
                      mtr.bootstraps = t(mtrEstimates),
                      te.ci1.90 = ci190,
                      te.ci1.95 = ci195,
                      te.ci2.90 = ci290,
                      te.ci2.95 = ci295,
                      mtr.ci1.90 = t(mtrci190),
                      mtr.ci1.95 = t(mtrci195),
                      mtr.ci2.90 = t(mtrci290),
                      mtr.ci2.95 = t(mtrci295),
                      prop.ci1.90 = t(propci190),
                      prop.ci1.95 = t(propci195),
                      prop.ci2.90 = t(propci290),
                      prop.ci2.95 = t(propci295))))

    } else if (point == FALSE) {
        stop("Bootstraps can only be applied if 'point == TRUE'.")
    } else {
        stop("Bootstraps must be an integer.")
    }
}


#' Estimation procedure from Mogstad, Torgovitsky (2017)
#'
#' This function estimates bounds on treatment effect parameters,
#' following the procedure described in Mogstad, Torgovitsky
#' (2017). Of the target parameters, the user can choose from the ATE,
#' ATT, ATU, LATE, and generalized LATE. The user is required to
#' provide a polynomial expression for the marginal treatment
#' responses (MTR), as well as a set of regressions. By restricting
#' the set of coefficients on each term of the MTRs to be consistent
#' with the regression estimates, the function is able to restrict
#' itself to a set of MTRs. The bounds on the treatment effect
#' parameter correspond to finding coefficients on the MTRs that
#' maximize their average difference.
#'
#' The estimation procedure relies on the propensity to take up
#' treatment. The propensity scores can either be estimated as part of
#' the estimation procedure, or the user can specify a variable in the
#' data set already containing the propensity scores.
#'
#' Constraints on the shape of the MTRs and marginal treatment effects
#' (MTE) can be imposed by the user, also. Specifically, bounds and
#' monotonicity restrictions are permitted. These constraints are only
#' enforced over a subset of the data. However, an audit procedure
#' randomly selects points outside of this subset to determine whether
#' or not the constraints hold. The user can specify how stringent
#' this audit procedure is using the function arguments.
#'
#' @param ivlike formula or vector of formulas used to specify the
#'     regressions for the IV-like estimands.
#' @param data \code{data.frame} used to estimate the treatment
#'     effects.
#' @param subset single subset condition or list of subset conditions
#'     corresponding to each IV-like estimand. See
#'     \code{\link{lists.mst}} on how to input the argument.
#' @param components a list of vectors of the terms/components from
#'     the regressions specifications we want to include in the set of
#'     IV-like estimands. To select the intercept term, include in the
#'     vector of variable names, `intercept'. See
#'     \code{\link{lists.mst}} on how to input the argument. If no
#'     argument is provided, then
#' @param propensity formula or variable name corresponding to
#'     propensity to take up treatment. If a formula is declared, then
#'     the function estimates propensity score according to the
#'     formula and link specified. If a variable name is declared,
#'     then the corresponding column in the data is taken as the
#'     vector of propensity scores.
#' @param link name of link function to estimate propensity score. Can
#'     be chosen from \code{linear}, \code{probit}, or \code{logit}.
#' @param treat variable name for treatment indicator
#' @param m0 one-sided formula for marginal treatment response
#'     function for control group. Splines can also be incorporated
#'     using the expression "uSplines(degree, knots, intercept)". The
#'     'intercept' argument may be omitted, and is set to \code{TRUE}
#'     by default.
#' @param m1 one-sided formula for marginal treatment response
#'     function for treated group. Splines can also be incorporated
#'     using the expression "uSplines(degree, knots, intercept)". The
#'     'intercept' argument may be omitted, and is set to \code{TRUE}
#'     by default.
#' @param uname variable name for unobservable used in declaring MTRs.
#' @param target target parameter to be estimated. Currently function
#'     allows for ATE ("\code{ate}"), ATT ("\code{att}"), ATU
#'     ("\code{atu}"), LATE ("\code{late}"), and generalized LATE
#'     ("\code{genlate}").
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
#' @param late.Z vector of variable names used to define the LATE.
#' @param late.from baseline set of values of Z used to define the
#'     LATE.
#' @param late.to comparison set of values of Z used to define the
#'     LATE.
#' @param late.X vector of variable names of covariates which we
#'     condition on when defining the LATE.
#' @param eval.X numeric vector of the values at which we condition
#'     variables in \code{late.X} on when estimating the LATE.
#' @param genlate.lb lower bound value of unobservable u for
#'     estimating generalized LATE.
#' @param genlate.ub upper bound value of unobservable u for
#'     estimating generalized LATE.
#' @param obseq.tol threshold for violation of observational
#'     equivalence. The threshold enters in multiplicatively. Thus, a
#'     value of 0 corresponds to no violation of observational
#'     equivalence, and the assumption that the model is correctly
#'     specified.
#' @param grid.nu number of evenly spread points in the interval [0,
#'     1] of the unobservable u used to form the grid for imposing
#'     shape restrictions on the MTRs.
#' @param grid.nx number of evenly spread points of the covariates to
#'     use to form the grid for imposing shape restrictions on the
#'     MTRs.
#' @param audit.nx number of points on the covariates space to audit
#'     in each iteration of the audit procedure.
#' @param audit.nu number of points in the interval [0, 1],
#'     corresponding to the normalized value of the unobservable term,
#'     to audit in each iteration of the audit procedure.
#' @param audit.max maximum number of iterations in the audit
#'     procedure.
#' @param audit.tol tolerance for determining when to end the audit
#'     procedure.
#' @param m1.ub numeric value for upper bound on MTR for treated
#'     group. By default, this will be set to the largest value of the
#'     observed outcome in the estimation sample.
#' @param m0.ub numeric value for upper bound on MTR for control
#'     group. By default, this will be set to the largest value of the
#'     observed outcome in the estimation sample.
#' @param m1.lb numeric value for lower bound on MTR for treated
#'     group. By default, this will be set to the smallest value of
#'     the observed outcome in the estimation sample.
#' @param m0.lb numeric value for lower bound on MTR for control
#'     group. By default, this will be set to the smallest value of
#'     the observed outcome in the estimation sample.
#' @param mte.ub numeric value for upper bound on treatment effect
#'     paramter of interest.
#' @param mte.lb numeric value for lower bound on treatment effect
#'     paramter of interest.
#' @param m0.dec logical, equal to TRUE if we want MTR for control
#'     group to be weakly monotone decreasing.
#' @param m0.inc logical, equal to TRUE if we want MTR for control
#'     group to be weakly monotone increasing.
#' @param m1.dec logical, equal to TRUE if we want MTR for treated
#'     group to be weakly monotone decreasing.
#' @param m1.inc logical, equal to TRUE if we want MTR for treated
#'     group to be weakly monotone increasing.
#' @param mte.dec logical, equal to TRUE if we want the MTE to be
#'     weakly monotone decreasing.
#' @param mte.inc logical, equal to TRUE if we want the MTE to be
#'     weakly monotone decreasing.
#' @param lpsolver name of the linear programming package in R used to
#'     obtain the bounds on the treatment effect.
#' @param point boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} if it is believed that the treatment effects are
#'     point identified. If set to \code{TRUE}, then a GMM procedure
#'     is implemented to estimate the treatment effects. Shape
#'     constraints on the MTRs will be ignored under point
#'     identification.
#' @param point.tol scalar, set default at 1-e08. Tolerance for bounds
#'     before automatically switchingto case of point
#'     identification. So if the estimated bounds are narrower than
#'     \code{point.tol}, and no shape restrictions are declared, the
#'     function instead provides a point estimate of the treatment
#'     effect. The output would be the same as if \code{point} was set
#'     to \code{TRUE}.
#' @param point.itermax integer, default of 2. Maximum number of
#'     iterations allowed for GMM estimation under point
#'     identification. So default estimate is the two-step GMM.
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @return Returns a list of results from throughout the estimation
#'     procedure. This includes all IV-like estimands; the propensity
#'     score model; bounds on the treatment effect; the estimated
#'     expectations of each term in the MTRs; the components and
#'     results of the LP problem.
#'
#' @examples
#' ivlikespecs <- c(ey ~ d | z,
#'                  ey ~ d | factor(z),
#'                  ey ~ d,
#'                  ey ~ d | factor(z))
#' jvec <- lists.mst(d, d, d, d)
#' svec <- lists.mst(, , , z %in% c(2, 4))
#'
#' mst(ivlike = ivlikespecs,
#'     data = dtm,
#'     components = jvec,
#'     propensity = pz,
#'     subset = svec,
#'     m0 = ~  u + u ^ 2,
#'     m1 = ~  u + u ^ 2,
#'     uname = u,
#'     target = "att",
#'     m0.dec = TRUE,
#'     m1.dec = TRUE)
#'
#' @export
ivmte.estimate <- function(ivlike, data, subset, components, propensity, link,
                  treat, m0, m1, uname = u, target, target.weight0,
                  target.weight1, target.knots0 = NULL, target.knots1 = NULL,
                  late.Z, late.from, late.to, late.X,
                  eval.X, genlate.lb, genlate.ub, obseq.tol = 1,
                  grid.nu = 20, grid.nx = 20, audit.nx = 2,
                  audit.nu = 3, audit.max = 5, audit.tol = 1e-08, m1.ub,
                  m0.ub, m1.lb, m0.lb, mte.ub, mte.lb, m0.dec, m0.inc,
                  m1.dec, m1.inc, mte.dec, mte.inc, lpsolver = NULL,
                  point = FALSE, point.itermax = 2,
                  point.tol = 1e-08, noisy = TRUE) {
    ## Match call arguments
    call <- match.call(expand.dots = FALSE)

    ## FIX: at some point, you may want to include "weights"

    ##---------------------------
    ## 0.a Check linear programming dependencies
    ##---------------------------

    if (is.null(lpsolver)) {
        if (requireNamespace("gurobi", quietly = TRUE)) {
            lpsolver <- "gurobi"
        } else if (requireNamespace("Rcplex", quietly = TRUE)) {
            lpsolver <- "Rcplex"
        } else if (requireNamespace("cplexAPI", quietly = TRUE)) {
            lpsolver <- "cplexAPI"
        } else if (requireNamespace("lpSolve", quietly = TRUE)) {
            lpsolver <- "lpSolve"
        } else {
            stop(gsub("\\s+", " ",
                      "Please install one of the following packages required for
                      estimation:
                      gurobi (version 7.5-1 or later);
                      Rcplex (version 0.3.3 or later);
                      cplexAPI (version 1.3.3 or later);
                      lpSolve (version 5.6.13 or later)."))
        }
    } else {
        if (! lpsolver %in% c("gurobi",
                              "Rcplex",
                              "cplexAPI",
                              "lpSolve",
                              "lpSolveAPI")) {
            stop(gsub("\\s+", " ",
                      paste0("Estimator is incompatible with linear programming
                             package '", lpsolver, "'. Please install one of the
                             following linear programming packages instead:
                             gurobi (version 7.5-1 or later);
                             Rcplex (version 0.3.3 or later);
                             cplexAPI (version 1.3.3 or later);
                             lpSolve (version 5.6.13 or later).")))
        }
    }

    ##---------------------------
    ## 0.b Check format of `formula', `subset', and `component' inputs
    ##---------------------------

    if (class_list(ivlike)) {

        ## Convert formula, components, and subset inputs into lists
        length_formula <- length(ivlike)

        userComponents  <- FALSE
        if (hasArg(components)) {
            if (!is.null(components)) {
                userComponents <- TRUE
            }
        }

        ## if (hasArg(components) & !is.null(components)) {
        if (userComponents) {
            ## userComponents <- TRUE
            if (class_list(components)) {
                length_components <- length(components)
                if (length_components == length_formula) {
                    specCompWarn <- TRUE
                } else {
                    specCompWarn <- FALSE
                }
            } else {
                length_components <- 1
                compList <- list()
                compList[[1]] <- components
                components <- complist
            }
        } else {
            ## userComponents <- FALSE
            length_components <- length_formula
            components <- as.list(replicate(length_formula, ""))
            warning(gsub("\\s+", " ",
                         "No list of components provided. All covariates in each
                         IV-like specification will be included when
                         constructing each S-set."),
                    call. = FALSE)
        }

        if (length_formula > length_components & length_components > 0) {
            warning(gsub("\\s+", " ",
                         "List of components not the same length of list of
                         IV-like specifications: more specifications than
                         component vectors. Specifications without corresponding
                         component vectors will include all covariates when
                         constructing the S-set."),
                    call. = FALSE)
            components[(length(components) + 1) : length(ivlike)] <- ""
        }

        if (length_formula < length_components) {
            warning(gsub("\\s+", " ",
                         "List of components not the same length of list of
                         IV-like specifications: more component vectors than
                         specifications. Component vectors without corresponding
                         specifications will be dropped."),
                    call. = FALSE)
            components <- components[1 : length(ivlike)]
        }

        ## Check the subset input---of the three lists that are input,
        ## only this can be omitted by the user, in which case no
        ## subsetting is used
        if (hasArg(subset)) {
            if (!class(subset) == "list") subset <- list(subset)

            ## Check if all subseting conditions are logical
            nonLogicSubset <- NULL
            for (i in 1:length(ivlike)) {
                if (subset[[i]] == "") {
                    ssubset <- replicate(nrow(data), TRUE)
                } else {
                    ssubset <- subset[[i]]
                }

                if (!is.logical(head(eval(substitute(ssubset), data)))) {
                    nonLogicSubset <- c(nonLogicSubset, i)
                }
            }
            if (length(nonLogicSubset) > 0) {
                stop(gsub("\\s+", " ",
                          paste0("The conditions in the following
                      positions of the subset list are not
                      logical: ",
                      paste(nonLogicSubset, collapse = ", "),
                      ". Please change the conditions so they
                      are logical.")))
            }

            if(length(subset) < length_formula) {
                warning(gsub("\\s+", " ",
                             "List of subset conditions not the same length
                              of list IV-like specifications: more
                              specifications than subsetting conditions.
                              Specifications without corresponding subset
                              conditions will include all observations."),
                        call. = FALSE)
                subset[length(subset) + 1 : length(ivlike)] <- ""
            }
            if(length(subset) > length_formula) {
                warning(gsub("\\s+", " ",
                             "List of subset conditions not the same length
                              of list IV-like specifications: more subset
                              conditions than IV-like specifications.
                              Subset conditions without corresponding
                              specifications will be dropped."),
                        call. = FALSE)
                subset <- subset[1 : length(ivlike)]
            }
        } else {
            ## if no subset input, then we construct it
            subset <- as.list(replicate(length_formula, ""))
        }

    } else {
        if (hasArg(subset)) {
            if (!is.logical(head(eval(substitute(subset), data)))) {
                stop(gsub("\\s+", " ",
                          "The subset condition is not logical.
                      Please change the condition to be logical."))
            }
        }

        specCompWarn <- FALSE
        if (hasArg(components)) {
            userComponents <- TRUE
            compString <- (unlist(lapply(components, deparse)))
            components <- components[compString != ""]
        } else {
            userComponents <- FALSE
            components <- as.list("")
            warning(gsub("\\s+", " ",
                         "No list of components provided. All covariates in each
                         IV-like specification will be included when
                         constructing each S-set."),
                    call. = FALSE)
        }
    }

    ##---------------------------
    ## 0.c Check numeric arguments and case completion
    ##---------------------------

    ## Variable and weight checks
    if (hasArg(treat)) {
        if (! deparse(substitute(treat)) %in% colnames(data)) {
            stop("Declared treatment indicator not found in data")
        }
    }

    if (hasArg(target)) {
        if (! target %in% c("ate", "att", "atu", "late", "genlate")) {
            stop(gsub("\\s+", " ",
                      "Specified target parameter is not recognized.
                      Choose from 'ate', 'att', 'atu', 'late', or 'genlate'."))
        }
        if (target == "late") {
            if (!(hasArg(late.Z) & hasArg(late.to) & hasArg(late.from))) {
                stop(gsub("\\s+", " ",
                          "Target paramter of 'late' requires arguments
                          'late.Z', 'late.to', and 'late.from'."))
            }

            if((hasArg(late.X) & !hasArg(eval.X)) |
               !hasArg(late.X) & hasArg(eval.X)) {
                stop(gsub("\\s+", " ",
                          "If the target parameter is 'late', then either both
                          late.X and eval.X are specified, or neither are
                          specified."))
            }
        }
        if (target == "genlate") {
            if (genlate.lb < 0 | genlate.ub > 1) {
                stop(gsub("\\s+", " ",
                          "'genlate.lb' and 'genlate.ub' must be between 0
                          and 1."))
            }
            if (genlate.lb >= genlate.ub) {
                stop("'genlate.lb' must be less than 'genlate.ub'.")
            }
        }
        if (target != "genlate" &
            (hasArg("genlate.lb") | hasArg("genlate.ub"))) {
            warning(gsub("\\s+", " ",
                         "Unless target parameter is 'genlate', 'genlate.lb' and
                         'genlate.ub' arguments will not be used."))
        }
        if (target != "late" & (hasArg("eval.X") | hasArg("late.X"))) {
            warning(gsub("\\s+", " ",
                         "Unless target parameter is 'late', 'eval.X' and
                         'late.X' arguments will not be used."))
        }
        if ((hasArg(target.weight0) | hasArg(target.weight1))) {
            stop(gsub("\\s+", " ",
                      "A preset target weight is chosen, and a custom target
                      weight is provided. Please provide an input for 'target'
                      only, or for 'target.weight0' and 'target.weight1'."))
        }
    } else {
        if (!(hasArg(target.weight0) & hasArg(target.weight1))) {
            stop(gsub("\\s+", " ",
                      "Only one target weight function is provided. If a custom
                      target weight is to be used, inputs for both
                      target.weight0 and target.weight1 must be provided."))
        }

        if (!hasArg(target)) {

            ## Convert all entries into lists
            target.weight0 <- c(target.weight0)
            target.weight1 <- c(target.weight1)
            target.knots0  <- c(target.knots0)
            target.knots1  <- c(target.knots1)

            ## Weight check
            weightCheck1 <- unlist(lapply(target.weight0, is.function))
            weightCheck2 <- unlist(lapply(target.weight0, function(x)
                is.numeric(x) * (length(x) == 1)))
            fails0 <- which(weightCheck1 + weightCheck2 == 0)

            if (length(fails0) > 0) {
                stop(gsub("\\s+", " ",
                          paste0("Each component of the custom weight
                             vectors/lists
                             must either be functions or constants. The
                             following entries of target.weight0 are neither: ",
                             paste(fails0, collapse = ", "),
                             ".")))
            }

            weightCheck1 <- unlist(lapply(target.weight1, is.function))
            weightCheck2 <- unlist(lapply(target.weight1, function(x)
                is.numeric(x) * (length(x) == 1)))
            fails1 <- which(weightCheck1 + weightCheck2 == 0)

            if (length(fails1) > 0) {
                stop(gsub("\\s+", " ",
                          paste0("Each component of the custom weight
                             vectors/lists
                             must either be functions or constants. The
                             following entries of target.weight1 are neither: ",
                             paste(fails1, collapse = ", "),
                             ".")))
            }

            if (length(target.weight0) != (length(target.knots0) + 1)) {
                stop(gsub("\\s+", " ",
                          paste0("The number of weight functions declared in
                                 target.weight0 must be exactly equal to the
                                 number of knots declared in target.knots0
                                 plus 1. Currently, the number of weights
                                 declared is ", length(target.weight0),
                                 ", and the number of knots declared is ",
                                 length(target.knots0), ".")))
            }
            if (length(target.weight1) != (length(target.knots1) + 1)) {
                stop(gsub("\\s+", " ",
                          paste0("The number of weight functions declared in
                                 target.weight1 must be exactly equal to the
                                 number of knots declared in target.knots1
                                 plus 1. Currently, the number of weights
                                 declared is ", length(target.weight1),
                                 ", and the number of knots declared is ",
                                 length(target.knots1), ".")))
            }

            ## Knots check
            if (!is.null(target.knots0)) {
                knotsCheck1 <- unlist(lapply(target.knots0, is.function))
                knotsCheck2 <- unlist(lapply(target.knots0, function(x)
                    is.numeric(x) * (length(x) == 1)))
                fails0 <- which(knotsCheck1 + knotsCheck2 == 0)

                if (length(fails0) > 0) {
                    stop(gsub("\\s+", " ",
                              paste0("Each component of the custom knots
                             vectors/lists
                             must either be functions or constants. The
                             following entries of target.knots0 are neither: ",
                             paste(fails0, collapse = ", "),
                             ".")))
                }
            }

            if (!is.null(target.knots1)) {
                knotsCheck1 <- unlist(lapply(target.knots1, is.function))
                knotsCheck2 <- unlist(lapply(target.knots1, function(x)
                    is.numeric(x) * (length(x) == 1)))
                fails1 <- which(knotsCheck1 + knotsCheck2 == 0)

                if (length(fails1) > 0) {
                    stop(gsub("\\s+", " ",
                              paste0("Each component of the custom knots
                             vectors/lists
                             must either be functions or constants. The
                             following entries of target.knots1 are neither: ",
                             paste(fails0, collapse = ", "),
                             ".")))
                }
            }
        }
    }

    if (hasArg(link)) {
        if (! link %in% c("linear", "logit", "probit")) {
            stop(gsub("\\s+", " ",
                      "Specified link is not recognized. Choose from 'linear',
                      'logit', or 'probit'."))
        }
    }

    if (hasArg(point)) {
        if (point == TRUE) {
            if (hasArg(m0.dec) | hasArg(m0.inc) |
                hasArg(m1.dec) | hasArg(m1.inc) |
                hasArg(mte.dec) | hasArg(mte.inc) |
                hasArg(audit.nu) | hasArg(audit.nx) |
                hasArg(grid.nu) | hasArg(grid.nx)|
                hasArg(audit.tol) | hasArg(audit.max)) {
                warning(gsub("\\s+", " ",
                             "If argument 'point' is set to TRUE, then shape
                             restrictions on m0 and m1 are ignored, and the
                             audit procedure is not implemented."))
            }
        }
    }

    ## Audit checks
    if (!(is.numeric(obseq.tol) & obseq.tol >= 0)) {
        stop("Cannot set obseq.tol below 0.")
    }
    if (!((grid.nu %% 1 == 0) & grid.nu >= 2)) {
        stop("grid.nu must be an integer greater than or equal to 2.")
    }
    if (!((grid.nx %% 1 == 0) & grid.nx >= 0)) {
        stop("grid.nx must be an integer greater than or equal to 0.")
    }

    if (!((audit.nx %% 1 == 0) & audit.nx > 0)) {
        stop("audit.nx must be an integer greater than or equal to 1.")
    }


    if (hasArg(m0.dec) | hasArg(m0.inc) |
        hasArg(m1.dec) | hasArg(m1.inc) |
        hasArg(mte.dec) | hasArg(mte.inc) |
        hasArg(m0.lb) | hasArg(m0.ub) |
        hasArg(m1.lb) | hasArg(m1.ub) |
        hasArg(mte.lb) | hasArg(mte.ub)) {

        noshape = FALSE ## indicator for whether shape restrictions declared

        if (!((audit.nu %% 1 == 0) & audit.nu > 0) | audit.nu < 2) {
            stop("audit.nu must be an integer greater than or equal to 2.")
        }
        
    } else {

        noshape = TRUE

        if (!((audit.nu %% 1 == 0) & audit.nu > 0)) {
            stop("audit.nu must be an integer greater than or equal to 1.")
        }
    }

    if (!((audit.max %% 1 == 0) & audit.max > 0)) {
        stop("audit.max must be an integer greater than or equal to 1.")
    }

    ##---------------------------
    ## 1. Restrict data to complete observations
    ##---------------------------

    ## Restrict data used for all parts of procedure to be the same.
    ## Collect list of all terms used in formula

    vars_y          <- NULL
    vars_formulas_x <- c()
    vars_formulas_z <- c()
    vars_subsets    <- c()
    vars_mtr        <- c()
    vars_weights    <- c()
    vars_propensity <- c()

    terms_formulas_x <- c()
    terms_formulas_z <- c()
    terms_mtr0       <- c()
    terms_mtr1       <- c()

    if (class_formula(ivlike)) {
        vars_formulas_x <- get_xz(ivlike)
        vars_formulas_z <- get_xz(ivlike, inst = TRUE)
        vars_y <- all.vars(ivlike)[1]
        terms_formulas <- attr(terms(Formula::as.Formula(ivlike)),
                               "term.labels")
    } else if (class_list(ivlike)) {
        if(!min(unlist(lapply(ivlike, class_formula)))) {
            stop(gsub("\\s+", " ",
                      "Not all elements in list of formulas are specified
                      correctly."))
        } else {
            vars_formulas_x <- unlist(lapply(ivlike,
                                             get_xz,
                                             inst = FALSE))
            vars_formulas_z <- unlist(lapply(ivlike,
                                             get_xz,
                                             inst = TRUE))

            vars_y <- unique(unlist(lapply(ivlike,
                                           function(x) all.vars(x)[[1]])))

            terms_formulas_x <- lapply(ivlike,
                                       get_xz,
                                       inst = FALSE,
                                       terms = TRUE)

            terms_formulas_z <- lapply(ivlike,
                                       get_xz,
                                       inst = TRUE,
                                       terms = TRUE)

            if (length(vars_y) > 1) {
                stop(gsub("\\s+", " ",
                          "Multiple response variables specified in list of
                          IV-like specifications."))
            }
        }
    } else {
        stop(gsub("\\s+", " ",
                  "'ivlike' argument must either be a formula or a vector of
                  formulas."))
    }

    ## Collect list of all terms in subsetting condition
    if (hasArg(subset)) {

        vnames <- colnames(data)

        if (class_list(subset)) {
            svec <- paste(unlist(lapply(subset, deparse)), collapse = " ")
        } else {
            svec <- deparse(substitute(subset))
        }
        svec <- subsetclean(svec)

        checklist <- c()
        checklist <- c(checklist, svec[!isfunctionstring(svec)])

        for (w in svec[isfunctionstring(svec)]) {
            arguments <- argstring(w)
            checklist <- c(checklist, unlist(strsplit(arguments, split = ",")))
        }

        for (w in checklist) {
            if (w %in% vnames) vars_subsets <- c(vars_subsets, w)
        }
    }

    ## Collect list of all terms used in MTRs
    if (class_formula(m1) & class_formula(m0)) {
        if(length(Formula::as.Formula(m0))[1] != 0 |
           length(Formula::as.Formula(m1))[1] != 0) {
            stop("m0 and m1 must be one-sided formulas.")
        }

        splinesobj <- list(removeSplines(m0),
                           removeSplines(m1))

        m0 <- splinesobj[[1]]$formula
        m1 <- splinesobj[[2]]$formula

        vars_mtr <- c(all.vars(splinesobj[[1]]$formula),
                      all.vars(splinesobj[[2]]$formula))

        if (!is.null(splinesobj[[1]]$formula)) {
            terms_mtr0 <- attr(terms(splinesobj[[1]]$formula),
                               "term.labels")
        }
        if (!is.null(splinesobj[[2]]$formula)) {
            terms_mtr1 <- attr(terms(splinesobj[[2]]$formula),
                               "term.labels")
        }

        if (!is.null(splinesobj[[1]]$splineslist)) {
            sf0 <- as.formula(paste("~",
                                    paste(unlist(splinesobj[[1]]$splineslist),
                                          collapse = " + ")))
            vars_mtr   <- c(vars_mtr, all.vars(sf0))
            terms_mtr0 <- c(terms_mtr0, attr(terms(sf0), "term.labels"))
        }

        if (!is.null(splinesobj[[2]]$splineslist)) {
            sf1 <- as.formula(paste("~",
                                    paste(unlist(splinesobj[[2]]$splineslist),
                                          collapse = " + ")))
            vars_mtr   <- c(vars_mtr, all.vars(sf1))
            terms_mtr1 <- c(terms_mtr1, attr(terms(sf1), "term.labels"))
        }

    } else {
        stop("m0 and m1 must be one-sided formulas.")
    }

    ## Collect list of variables used in custom weights (if defined)
    if (!hasArg(target)) {

        for (i in 1:length(target.weight0)) {
            if (is.function(target.weight0[[i]])) {
                vars_weights <- c(vars_weights,
                                  formalArgs(target.weight0[[i]]))
            }

            if (i < length(target.weight0)) {
                if (is.function(target.knots0[[i]])) {
                    vars_weights <- c(vars_weights,
                                      formalArgs(target.knots0[[i]]))
                }
            }
        }

        for (i in 1:length(target.weight1)) {
            if (is.function(target.weight1[[i]])) {
                vars_weights <- c(vars_weights,
                                  formalArgs(target.weight1[[i]]))
            }

            if (i < length(target.weight1)) {
                if (is.function(target.knots1[[i]])) {
                    vars_weights <- c(vars_weights,
                                      formalArgs(target.knots1[[i]]))
                }
            }
        }
    }

    ## Collect list of all terms used in propensity formula

    if (hasArg(propensity)) {
        if (class_formula(propensity)) {
            ptreat <- all.vars(propensity)[1]
            vars_propensity <- all.vars(propensity)

            if (hasArg(treat)) {
                if (ptreat != deparse(substitute(treat))) {
                    warning(gsub("\\s+", " ",
                                 "Variable listed in 'treat' argument differs
                                 from dependent variable in propensity score
                                 formula. Dependent variable from propensity
                                 score formula will be used as the treatment
                                 variable."),
                            call. = FALSE)
                    treat <- ptreat
                } else {
                    treat <- deparse(substitute(treat))
                }
            } else {
                warning(gsub("\\s+", " ",
                             paste0("'treat' argument is not declared.
                              Dependent variable from the propensity
                              score formula, '",
                              ptreat,
                              "', will be used as the treatment variable.")))
                treat <- ptreat
            }

        } else {

            if (! deparse(substitute(propensity)) %in% colnames(data)) {
                stop(gsub("\\s+", " ",
                          "Propensity score argument is interpreted as a
                          variable name, but is not found in the data set."))
            }

            vars_propensity <- c(vars_propensity,
                                 deparse(substitute(propensity)))

            ## Determine treatment variable
            if (hasArg(treat)) {
                treat <- deparse(substitute(treat))
                vars_propensity <- c(vars_propensity,
                                      treat)
            } else if (class(ivlike) == "formula") {
                warning(gsub("\\s+", " ",
                             "First independent variable of IV-like
                             specification regression is selected as the
                             treatment variable."),
                        call. = FALSE)
                treat <- all.vars(ivlike)[2]
                vars_propensity <- c(vars_propensity,
                                      treat)
            } else if (is.list(ivlike)) {
                warning(gsub("\\s+", " ",
                             "First independent variable of first IV regression
                             is selected as the treatment variable."),
                        call. = FALSE)
                treat <- all.vars(ivlike[[1]])[2]
                vars_propensity <- c(vars_propensity,
                                      treat)
            } else {
                stop("Treatment variable indeterminable.")
            }
        }
    } else {
        ## Determine treatment variable
        if (hasArg(treat)) {
            treat <- deparse(substitute(treat))
            vars_propensity <- treat
        } else if (class(ivlike) == "formula") {
            warning(gsub("\\s+", " ",
                         "First independent variable of IV-like
                             specification regression is selected as the
                             treatment variable."))
            treat <- all.vars(ivlike)[2]
        } else if (is.list(ivlike)) {
            warning(gsub("\\s+", " ",
                         "First independent variable of first IV regression
                             is selected as the treatment variable."))
            treat <- all.vars(ivlike[[1]])[2]
        } else {
            stop("Treatment variable indeterminable.")
        }

        ## Construct propensity formula
        terms_propensity <- c(unlist(terms_formulas_x),
                              unlist(terms_formulas_z),
                              unlist(terms_mtr0),
                              unlist(terms_mtr1))

        ## Remove all u terms
        uterms <- c()

        um1 <- which(terms_propensity == deparse(substitute(uname)))
        um2 <- grep(paste0("^[", deparse(substitute(uname)), "][[:punct:]]"),
                    terms_propensity)
        um3 <- grep(paste0("[[:punct:]][", deparse(substitute(uname)), "]$"),
                    terms_propensity)
        um4 <- grep(paste0("[[:punct:]][",
                           deparse(substitute(uname)),
                           "][[:punct:]]"),
                    terms_propensity)
        um5 <- grep(paste0("[[:punct:]][", deparse(substitute(uname)), "]\\s+"),
                    terms_propensity)
        um6 <- grep(paste0("\\s+[", deparse(substitute(uname)), "][[:punct:]]"),
                    terms_propensity)
        um7 <- grep(paste0("\\s+[", deparse(substitute(uname)), "]\\s+"),
                    terms_propensity)
        um8 <- grep("uSplines", terms_propensity)

        for (i in 1:8) {
            uterms <- c(uterms, terms_propensity[get(paste0("um", i))])
        }
        terms_propensity <- terms_propensity[!(terms_propensity %in% uterms) &
                                           !(terms_propensity == treat)]
        propensity <- paste(treat,
                            paste(unique(terms_propensity),
                                  collapse = " + "),
                            sep = " ~ ")
        propWarn <- paste("No propensity score formula or propensity score
                           variable name provided. By default, the function will
                           create a propensity score formula using all the
                           covariates in the IV-like specifications and the
                           MTRs. The propensity score formula generated is:",
                           propensity)
        warning(gsub("\\s+", " ", propWarn), call. = FALSE)
        propensity <- as.formula(propensity)
    }

    ## Remove unobserved variable from list
    allvars <- unique(c(vars_y,
                        vars_formulas_x,
                        vars_formulas_z,
                        vars_subsets,
                        vars_mtr,
                        vars_weights,
                        vars_propensity))
    allvars <- allvars[allvars != deparse(substitute(uname))]

    newpropensity <- unique(c(vars_formulas_x,
                              vars_formulas_z,
                              vars_mtr))

    newpropensity <- newpropensity[(newpropensity !=
                                    deparse(substitute(uname))) &
                                   (newpropensity != treat)]

    comp_filler <- lapply(terms_formulas_x,
                          function(x) as.character(unstring(x)))

    ## Fill in components list if necessary
    if (userComponents) {
        compMissing1 <- unlist(lapply(components, function(x) deparse(x) == ""))
        compMissing2 <- unlist(lapply(components, function(x) x == ""))
        compMissing <- as.logical(compMissing1 + compMissing2)
        if (sum(compMissing) > 0 & specCompWarn) {
            warning(gsub("\\s+", " ",
                         "Specifications without corresponding
                         component vectors will include all covariates when
                         constructing the S-set."),
                    call. = FALSE)
        }
        components[compMissing] <- comp_filler[compMissing]
    } else {
        components <- comp_filler
    }

    ## Keep only complete cases

    varError <- allvars[! allvars[allvars != "intercept"] %in%
                        colnames(data)]
    if (length(varError) > 0) {
        varError <- paste0("The following variables are not contained
                          in the data set: ",
                          paste(varError, collapse = ", "),
                          ".")
        stop(gsub("\\s+", " ", varError), call. = FALSE)
    }

    data  <- data[(complete.cases(data[, allvars[allvars != "intercept"]])), ]
    ## Adjust row names to handle bootstrapping
    rownames(data) <- as.character(seq(1, nrow(data)))
    cdata <- data

    ##---------------------------
    ## 2. Obtain propensity scores
    ##---------------------------

    if (noisy == TRUE) {
        message("Obtaining propensity scores...\n")
    }

    ## Estimate propensity scores
    if (class_formula(propensity)) {
        pcall <- modcall(call,
                         newcall = propensity.mst,
                         keepargs = c("link", "late.Z", "late.X"),
                         dropargs = "propensity.mst",
                         newargs = list(data = quote(cdata),
                                        formula = propensity))
    } else {
        pcall <- modcall(call,
                         newcall = propensity.mst,
                         keepargs = c("link", "late.Z", "late.X"),
                         dropargs = "propensity.mst",
                         newargs = list(data = quote(cdata),
                                        formula = substitute(propensity)))
    }
    pmodel <- eval(pcall)

    ##---------------------------
    ## 3. Generate target moments/gamma terms
    ##---------------------------

    xindex0 <- NULL
    xindex1 <- NULL
    uexporder0 <- NULL
    uexporder1 <- NULL

    if (noisy == TRUE) {
        message("Generating target moments...\n")
    }

    if (!is.null(m0)) {
        m0call <- modcall(call,
                          newcall = polyparse.mst,
                          keepargs = c("uname"),
                          newargs = list(formula = m0,
                                         data = quote(cdata)))
    } else {
        m0call <- NULL
    }

    if (!is.null(m1)) {
        m1call <- modcall(call,
                          newcall = polyparse.mst,
                          keepargs = c("uname"),
                          newargs = list(formula = m1,
                                         data = quote(cdata)))
    } else {
        m1call <- NULL
    }

    ## Generate target weights
    if (hasArg(target)) {
        if (target == "ate") {
            w1 <- wate1.mst(cdata)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "att") {
            w1 <- watt1.mst(cdata, mean(cdata[[treat]]), pmodel$phat)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "atu") {
            w1 <- watu1.mst(cdata, 1 - mean(cdata[[treat]]), pmodel$phat)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "late") {
            if (!hasArg(late.X)) {
                late.X <- NULL
                eval.X <- NULL
            }
            w1 <- wlate1.mst(cdata, late.from, late.to, substitute(late.Z),
                             pmodel$model, substitute(late.X), eval.X)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "genlate") {
            w1 <- wgenlate1.mst(cdata, genlate.lb, genlate.ub)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else {
            stop("Unrecognized target parameter.")
        }

        ## Integrate m0 and m1 functions
        if (!is.null(m0)) {
            if (noisy == TRUE) {
                message("    Integrating terms for control group...\n")
            }
            pm0 <- eval(as.call(m0call))

            if (point == FALSE) {
                gstar0 <- gengamma.mst(pm0, w0$lb, w0$ub, w0$mp)
            } else {
                gstar0 <- gengamma.mst(pm0, w0$lb, w0$ub, w0$mp, means = FALSE)
                xindex0 <- c(xindex0, pm0$xindex)
                uexporder0 <- c(uexporder0, pm0$exporder)
            }
        } else {
            gstar0 <- NULL
            pm0 <- NULL
        }

        if (!is.null(m1)) {
            if (noisy == TRUE) {
                message("    Integrating terms for treated group...\n")
            }
            pm1 <- eval(as.call(m1call))

            if (point == FALSE) {
                gstar1 <- gengamma.mst(pm1, w1$lb, w1$ub, w1$mp)
            } else {
                gstar1 <- gengamma.mst(pm1, w1$lb, w1$ub, w1$mp, means = FALSE)
                xindex1 <- c(xindex1, pm1$xindex)
                uexporder1 <- c(uexporder1, pm1$exporder)
            }
        } else {
            gstar1 <- NULL
            pm1 <- NULL
        }

        if (point == FALSE) {
            gstarSplineObj0 <- genGammaSplines.mst(splines = splinesobj[[1]],
                                                data = cdata,
                                                lb = w0$lb,
                                                ub = w0$ub,
                                                multiplier = w0$mp,
                                                d = 0)
            gstarSpline0 <- gstarSplineObj0$gamma

            gstarSplineObj1 <- genGammaSplines.mst(splines = splinesobj[[2]],
                                                data = cdata,
                                                lb = w1$lb,
                                                ub = w1$ub,
                                                multiplier = w1$mp,
                                                d = 1)
            gstarSpline1 <- gstarSplineObj1$gamma
        } else {
            gstarSplineObj0 <- genGammaSplines.mst(splines = splinesobj[[1]],
                                                data = cdata,
                                                lb = w0$lb,
                                                ub = w0$ub,
                                                multiplier = w0$mp,
                                                d = 0,
                                                means = FALSE)
            gstarSpline0 <- gstarSplineObj0$gamma
            xindex0 <- c(xindex0, gstarSplineObj0$interactions)
            uexporder0 <- c(uexporder0,
                            rep(-1, length(xindex0) - length(uexporder0)))

            gstarSplineObj1 <- genGammaSplines.mst(splines = splinesobj[[2]],
                                                data = cdata,
                                                lb = w1$lb,
                                                ub = w1$ub,
                                                multiplier = w1$mp,
                                                d = 1,
                                                means = FALSE)
            gstarSpline1 <- gstarSplineObj1$gamma
            xindex1 <- c(xindex1, gstarSplineObj1$interactions)
            uexporder1 <- c(uexporder1,
                            rep(-1, length(xindex1) - length(uexporder1)))
        }
    } else {

        ## Convert fixed/numeric weights into functions
        if (is.numeric(target.weight0)) {
            target.weight0 <- sapply(target.weight0, constructConstant)
        } else {
            numeric <- which(unlist(lapply(target.weight0, is.numeric)))
            target.weight0[numeric] <- sapply(unlist(target.weight0[numeric]),
                                              constructConstant)
        }

        if (is.numeric(target.weight1)) {
            target.weight1 <- sapply(target.weight1, constructConstant)
        } else {
            numeric <- which(unlist(lapply(target.weight1, is.numeric)))
            target.weight1[numeric] <- sapply(unlist(target.weight1[numeric]),
                                              constructConstant)
        }

        ## Convert fixed/numeric knots into functions
        if (!is.null(target.knots0)) {
            if (is.numeric(target.knots0)) {
                target.knots0 <- sapply(target.knots0, constructConstant)
            } else {
                numeric <- which(unlist(lapply(target.knots0, is.numeric)))
                target.knots0[numeric] <- sapply(unlist(target.knots0[numeric]),
                                                 constructConstant)
            }
        }

        if (!is.null(target.knots1)) {
            if (is.numeric(target.knots1)) {
                target.knots1 <- sapply(target.knots1, constructConstant)
            } else {
                numeric <- which(unlist(lapply(target.knots1, is.numeric)))
                target.knots1[numeric] <- sapply(unlist(target.knots1[numeric]),
                                                 constructConstant)
            }
        }

        for (d in 0:1) {

            mtr <- get(paste0("m", d))

            ## Include end points
            splitFirst <- function(...) {
                0
            }

            splitLast <- function(...) {
                1
            }

            target.knots <- c(splitFirst,
                              get(paste0("target.knots", d)),
                              splitLast)

            target.weight <- get(paste0("target.weight", d))

            ## Integrate non-splines terms

            if (!is.null(mtr)) {
                if (noisy == TRUE) {
                    if (d == 0) {
                        message(
                      "    Integrating non-spline terms for control group...\n")
                    } else {
                        message(
                      "    Integrating non-spline terms for treated group...\n")
                    }
                }

                pm <- eval(as.call(get(paste0("m", d, "call"))))

                gamma <- rep(0, length(pm$terms))

                for(i in 1:length(target.weight)) {

                    wKnotVarsL <- formalArgs(target.knots[[i]])
                    wKnotVarsU <- formalArgs(target.knots[[i + 1]])

                    if (wKnotVarsL[1] == "...") {
                        lb <- unlist(lapply(X = seq(1, nrow(cdata)),
                                            FUN = funEval,
                                            fun = target.knots[[i]]))
                    } else {
                        lb <- unlist(lapply(X = split(cdata[, wKnotVarsL],
                                                      seq(1, nrow(cdata))),
                                            FUN = funEval,
                                            fun = target.knots[[i]],
                                            argnames = wKnotVarsL))
                    }

                    if (wKnotVarsU[1] == "...") {
                        ub <- unlist(lapply(X = seq(1, nrow(cdata)),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]]))

                    } else {
                        ub <- unlist(lapply(X = split(cdata[, wKnotVarsU],
                                                      seq(1, nrow(cdata))),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]],
                                            argnames = wKnotVarsU))
                    }

                    wValVars  <- formalArgs(target.weight[[i]])

                    if (wValVars[1] == "...") {
                        weights <- unlist(lapply(X = seq(1, nrow(cdata)),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]]))
                    } else {
                        weights <- unlist(lapply(X = split(cdata[, wValVars],
                                                           seq(1, nrow(cdata))),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]],
                                                 argnames = wValVars))
                    }

                    if (point == FALSE) {
                        gamma <- gamma + gengamma.mst(pm,
                                                      lb = lb,
                                                      ub = ub,
                                                      multiplier = weights)
                    } else {
                        gamma <- gamma + gengamma.mst(pm,
                                                      lb = lb,
                                                      ub = ub,
                                                      multiplier = weights,
                                                      means = FALSE)
                    }
                }

                assign(paste0("gstar", d), gamma)
                assign(paste0("pm", d), pm)
            } else {
                assign(paste0("gstar", d), NULL)
                assign(paste0("pm", d), NULL)
            }

            ## Integrate splines terms
            if (noisy == TRUE) {
                if (d == 0) {
                    message(
                        "    Integrating spline terms for control group...\n")
                } else {
                    message(
                        "    Integrating spline terms for treated group...\n")
                }
            }

            noSplineMtr <- splinesobj[[d + 1]]

            if (!is.null(noSplineMtr$splineslist)) {

                gammaSplines <- 0

                for (i in 1:length(target.weight)) {

                    wKnotVarsL <- formalArgs(target.knots[[i]])
                    wKnotVarsU <- formalArgs(target.knots[[i + 1]])

                    if (wKnotVarsL[1] == "...") {
                        lb <- unlist(lapply(X = seq(1, nrow(cdata)),
                                            FUN = funEval,
                                            fun = target.knots[[i]]))
                    } else {
                        lb <- unlist(lapply(X = split(cdata[, wKnotVarsL],
                                                      seq(1, nrow(cdata))),
                                            FUN = funEval,
                                            fun = target.knots[[i]],
                                            argnames = wKnotVarsL))
                    }

                    if (wKnotVarsU[1] == "...") {
                        ub <- unlist(lapply(X = seq(1, nrow(cdata)),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]]))

                    } else {
                        ub <- unlist(lapply(X = split(cdata[, wKnotVarsU],
                                                      seq(1, nrow(cdata))),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]],
                                            argnames = wKnotVarsU))
                    }

                    wValVars  <- formalArgs(target.weight[[i]])

                    if (wValVars[1] == "...") {
                        weights <- unlist(lapply(X = seq(1, nrow(cdata)),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]]))
                    } else {
                        weights <- unlist(lapply(X = split(cdata[, wValVars],
                                                           seq(1, nrow(cdata))),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]],
                                                 argnames = wValVars))
                    }

                    if (point == FALSE) {
                        gammaSplines <- gammaSplines +
                            genGammaSplines.mst(splines = noSplineMtr,
                                                data = cdata,
                                                lb = lb,
                                                ub = ub,
                                                multiplier = weights,
                                                d = d)
                    } else {
                        gammaSplines <- gammaSplines +
                            genGammaSplines.mst(splines = noSplineMtr,
                                                data = cdata,
                                                lb = lb,
                                                ub = ub,
                                                multiplier = weights,
                                                d = d,
                                                means = FALSE)
                    }
                }
                assign(paste0("gstarSpline", d), gammaSplines)
            } else {
                assign(paste0("gstarSpline", d), NULL)
            }
        }
    }

    if (point == FALSE) {
        gstar0 <- c(gstar0, gstarSpline0)
        gstar1 <- c(gstar1, gstarSpline1)
    } else {
        gstar0 <- cbind(gstar0, gstarSpline0)
        gstar1 <- cbind(gstar1, gstarSpline1)
    }

    ##---------------------------
    ## 4. Generate moments/gamma terms for IV-like estimands
    ##---------------------------

    if (noisy == TRUE) {
        message("Generating IV-like moments...")
    }

    sset  <- list() ## Contains all IV-like estimates and their
                    ## corresponding moments/gammas
    scount <- 1     ## counter for S-set constraints

    ## Construct `sset' object when a single IV-like specification is
    ## provided
    if (class_formula(ivlike)) {

        ## Obtain coefficient estimates and S-weights
        scall <- modcall(call,
                       newcall = ivlike.mst,
                       keepargs = c("subset"),
                       newargs = list(formula = ivlike,
                                      treat = quote(treat),
                                      data = quote(cdata),
                                      components = components))

        sest <- eval(scall)

        ncomponents <- length(sest$betas)

        if (hasArg(subset)) {
            subset_index <- rownames(data[eval(substitute(subset), data), ])
        } else {
            subset_index <- rownames(data)
        }

        ## Generate moments (gammas) corresponding to IV-like
        ## estimands
        if (point == FALSE) {
            setobj <- gensset.mst(data = cdata,
                                  sset = sset,
                                  sest = sest,
                                  splinesobj = splinesobj,
                                  pmodobj = pmodel$phat[subset_index],
                                  pm0 = pm0,
                                  pm1 = pm1,
                                  ncomponents = ncomponents,
                                  scount = scount,
                                  subset_index = subset_index,
                                  noisy = noisy)
        } else {
            setobj <- gensset.mst(data = cdata,
                                  sset = sset,
                                  sest = sest,
                                  splinesobj = splinesobj,
                                  pmodobj = pmodel$phat[subset_index],
                                  pm0 = pm0,
                                  pm1 = pm1,
                                  ncomponents = ncomponents,
                                  scount = scount,
                                  subset_index = subset_index,
                                  means = FALSE,
                                  yvar = vars_y,
                                  dvar = treat,
                                  noisy = noisy)
        }

        sset <- setobj$sset
        scount <- setobj$scount

    } else if (class_list(ivlike)) {
        ## Construct `sset' object when multiple IV-like
        ## specifications are provided

        ## loop across IV specifications
        for (i in 1:length(ivlike)) {

            sformula   <- ivlike[[i]]
            scomponent <- components[[i]]

            if (subset[[i]] == "") {
                ssubset <- replicate(nrow(data), TRUE)
            } else {
                ssubset <- subset[[i]]
            }

            ## Obtain coefficient estimates and S-weights
            ## corresponding to the IV-like estimands
            sdata <- data[eval(substitute(ssubset), data), ]
            sest  <- ivlike.mst(formula = sformula,
                                data = sdata,
                                components = scomponent,
                                treat = treat,
                                list = TRUE)

            ## Generate moments (gammas) corresponding to IV-like
            ## estimands
            subset_index <- rownames(sdata)
            ncomponents <- length(sest$betas)
            pmodobj <- pmodel$phat[subset_index]
            if (point == FALSE) {
                setobj <- gensset.mst(data = cdata,
                                      sset = sset,
                                      sest = sest,
                                      splinesobj = splinesobj,
                                      pmodobj = pmodobj,
                                      pm0 = pm0,
                                      pm1 = pm1,
                                      ncomponents = ncomponents,
                                      scount = scount,
                                      subset_index = subset_index,
                                      noisy = noisy)
            } else {
                setobj <- gensset.mst(data = cdata,
                                      sset = sset,
                                      sest = sest,
                                      splinesobj = splinesobj,
                                      pmodobj = pmodobj,
                                      pm0 = pm0,
                                      pm1 = pm1,
                                      ncomponents = ncomponents,
                                      scount = scount,
                                      subset_index = subset_index,
                                      means = FALSE,
                                      yvar = vars_y,
                                      dvar = treat,
                                      noisy = noisy)
            }

            ## Update set of moments (gammas)
            sset <- setobj$sset
            scount <- setobj$scount
        }
    } else {
        stop(gsub("\\s+", " ",
                  "'ivlike' argument must either be a formula or a vector of
                  formulas."))
    }

    ## Prepare GMM estimate estimate if `point' agument is set to TRUE
    if (point == TRUE) {

        ## Obtain GMM estimate

        gmmResult <- gmmEstimate(sset = sset,
                                 gstar0 = gstar0,
                                 gstar1 = gstar1,
                                 itermax = point.itermax,
                                 tol = point.tol,
                                 noisy = noisy)

        return(list(sset  = sset,
                    gstar = list(g0 = gstar0,
                                 g1 = gstar1),
                    propensity = pmodel,
                    te = gmmResult$te,
                    mtr.coef = gmmResult$coef))
    }

    ##---------------------------
    ## 5. Define constraint matrices using the audit
    ##---------------------------

    if (noisy == TRUE) {
        if (obseq.tol > 0) {
            message("\nPerforming audit procedure...\n")
        }
    }

    audit.args <- c("uname", "grid.nu", "grid.nx",
                    "audit.nx", "audit.nu", "audit.max", "audit.tol",
                    "m1.ub", "m0.ub",
                    "m1.lb", "m0.lb",
                    "mte.ub", "mte.lb", "m0.dec",
                    "m0.inc", "m1.dec", "m1.inc", "mte.dec",
                    "mte.inc", "obseq.tol")

    audit_call <- modcall(call,
                          newcall = audit.mst,
                          keepargs = audit.args,
                          newargs = list(data = quote(cdata),
                                         m0   = quote(m0),
                                         m1   = quote(m1),
                                         splinesobj = quote(splinesobj),
                                         vars_mtr = quote(vars_mtr),
                                         terms_mtr0 = quote(terms_mtr0),
                                         terms_mtr1 = quote(terms_mtr1),
                                         sset = quote(sset),
                                         gstar0 = quote(gstar0),
                                         gstar1 = quote(gstar1),
                                         lpsolver = quote(lpsolver)))

    ## Impose default upper and lower bounds on m0 and m1
    if (!hasArg(m1.ub) | !hasArg(m0.ub)) {
        maxy <- max(cdata[, vars_y])
        if (!hasArg(m1.ub)) {
            audit_call <- modcall(audit_call,
                                  newargs = list(m1.ub = maxy,
                                                 m1.ub.default = TRUE))
        }
        if (!hasArg(m0.ub)) {
            audit_call <- modcall(audit_call,
                                  newargs = list(m0.ub = maxy,
                                                 m0.ub.default = TRUE))
        }
    }
    if (!hasArg(m1.lb) | !hasArg(m0.lb)) {
        miny <- min(cdata[, vars_y])
        if (!hasArg(m1.lb)) {
            audit_call <- modcall(audit_call,
                                  newargs = list(m1.lb = miny,
                                                 m1.lb.default = TRUE))
        }
        if (!hasArg(m0.lb)) {
            audit_call <- modcall(audit_call,
                                  newargs = list(m0.lb = miny,
                                                 m0.lb.default = TRUE))
        }
    }

    ##---------------------------
    ## 6. Obtain the bounds
    ##---------------------------

    audit <- eval(audit_call)

    ## stop("end of test")
    
    ## If there are no shape restrictions and bounds are very tight,
    ## apply point identification method
    if (abs(audit$max - audit$min) < point.tol & noshape == TRUE) {

        message(gsub("\\s+", " ",
                     paste0("Width of bounds is narrower than 'point.tol', which
                     is set to ", point.tol, ". Estimation under point
                     identification will instead be performed.")))


        stop("You need to reconstruct all the ssets and gstar objects if the bounds are so small that you switch to the the GMM estimate.")
        gmmResult <- gmmEstimate(sset = sset, gstar0 = gstar0, gstar1
        = gstar1, itermax = point.itermax, tol = point.tol, noisy =
        noisy)

        return(list(sset  = sset,
                    gstar = list(g0 = gstar0,
                                 g1 = gstar1),
                    propensity = pmodel,
                    te = gmmResult$te,
                    se = gmmResult$se,
                    ci90 = gmmResult$ci90,
                    ci95 = gmmResult$ci95,
                    mtr.coef = gmmResult$coef,
                    mtr.vcov = gmmResult$vcov))

    } else {
        cat("Bound: (", audit$min, ",", audit$max, ")\n")

        ## include additional output material
        return(list(sset  = sset,
                    gstar = list(g0 = gstar0,
                                 g1 = gstar1),
                    propensity = pmodel,
                    bound = c(audit$min, audit$max),
                    lpresult =  audit$lpresult,
                    poly0 = pm0,
                    poly1 = pm1,
                    auditgrid = audit$gridobj,
                    minobseq = audit$minobseq,
                    splinesdict = list(splinesobj[[1]]$splinesdict,
                                       splinesobj[[2]]$splinesdict)))
    }
}


#' Generating LP moments for IV-like estimands
#'
#' This function takes in the IV estimate and its IV-like
#' specification, and generates a list containing the corresponding
#' point estimate, and the corresponding moments (gammas) that will
#' enter into the constraint matrix of the LP problem.
#'
#' @param data \code{data.frame} used to estimate the treatment
#'     effects.
#' @param sset A list, which is modified and returned as the output.
#' @param sest A list containing the point estimates and S-weights
#'     corresponding to a particular IV-like estimand.
#' @param splinesobj list of spline components in the MTRs for treated
#'     and control groups. Spline terms are extracted using
#'     \code{\link{removeSplines}}.
#' @param pmodobj A vector of propensity scores.
#' @param pm0 A list of the monomials in the MTR for d = 0.
#' @param pm1 A list of the monomials in the MTR for d = 1.
#' @param ncomponents The number of components from the IV regression
#'     we want to include in the S-set.
#' @param scount A counter for the number of elements in the S-set.
#' @param subset_index An index for the subset of the data the IV
#'     regression is restricted to.
#' @param means boolean, set to \code{TRUE} by default. If set to
#'     \code{TRUE}, then the gamma moments are returned, i.e. sample
#'     averages are taken. If set to \code{FALSE}, then no sample
#'     averages are taken, and a matrix is returned. The sample
#'     average of each column of the matrix corresponds to a
#'     particular gamma moment.
#' @param yvar name of outcome variable. This is only used if
#'     \code{means = FALSE}, which occurs when the user believes the
#'     treatment effect is point identified.
#' @param dvar name of treatment indicator. This is only used if
#'     \code{means = FALSE}, which occurs when the user believes the
#'     treatment effect is point identified.
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @return A list containing the point estimate for the IV regression,
#'     and the expectation of each monomial term in the MTR.
gensset.mst <- function(data, sset, sest, splinesobj, pmodobj, pm0, pm1,
                        ncomponents, scount, subset_index, means = TRUE,
                        yvar, dvar, noisy = TRUE) {

    for (j in 1:ncomponents) {
        if (noisy == TRUE) {
            message(paste0("    Moment ", scount, "..."))
        }

        if (!is.null(pm0)) {
            if (means == TRUE) {
                gs0 <- gengamma.mst(monomials = pm0,
                                    lb = pmodobj,
                                    ub = 1,
                                    multiplier = sest$sw0[, j],
                                    subset = subset_index)
            } else {
                gs0 <- gengamma.mst(monomials = pm0,
                                    lb = pmodobj,
                                    ub = 1,
                                    multiplier = sest$sw0[, j],
                                    subset = subset_index,
                                    means = FALSE)

                ## print(paste("S-weight of", j, "th componnt from S set, d = 0"))
                ## print(sest$sw0[, j])
            }
        } else {
            gs0 <- NULL
        }

        if (!is.null(pm1)) {
            if (means == TRUE) {
                gs1 <- gengamma.mst(monomials = pm1,
                                    lb = 0,
                                    ub = pmodobj,
                                    multiplier = sest$sw1[, j],
                                    subset = subset_index)
            } else {
                gs1 <- gengamma.mst(monomials = pm1,
                                    lb = 0,
                                    ub = pmodobj,
                                    multiplier = sest$sw1[, j],
                                    subset = subset_index,
                                    means = FALSE)
                ## print(paste(j, "th componnt from S set, d = 1"))
                ## print(sest$sw1[, j])
            }

            ## print("S-weight ratios, (d = 0) / (d = 1)")
            ## print(sest$sw0[, j] / sest$sw1[, j])
            ## print(paste(j, "th componnt from S set, d = 1"))
            ## print(head(sest$sw1[, j]))

        } else {
            gs1 <- NULL
        }
        
        if (means == TRUE) {
            gsSpline0 <- genGammaSplines.mst(splines = splinesobj[[1]],
                                             data = data,
                                             lb = pmodobj,
                                             ub = 1,
                                             multiplier = sest$sw0[, j],
                                             subset = subset_index,
                                             d = 0)$gamma

            gsSpline1 <- genGammaSplines.mst(splines = splinesobj[[2]],
                                             data = data,
                                             lb = 0,
                                             ub = pmodobj,
                                             multiplier = sest$sw1[, j],
                                             subset = subset_index,
                                             d = 1)$gamma
        } else {
            gsSpline0 <- genGammaSplines.mst(splines = splinesobj[[1]],
                                             data = data,
                                             lb = pmodobj,
                                             ub = 1,
                                             multiplier = sest$sw0[, j],
                                             subset = subset_index,
                                             d = 0,
                                             means = FALSE)$gamma

            gsSpline1 <- genGammaSplines.mst(splines = splinesobj[[2]],
                                             data = data,
                                             lb = 0,
                                             ub = pmodobj,
                                             multiplier = sest$sw1[, j],
                                             subset = subset_index,
                                             d = 1,
                                             means = FALSE)$gamma
        }

        ## generate components of constraints
        if (means == TRUE) {
            sset[[paste0("s", scount)]] <- list(beta = sest$beta[j],
                                                g0 = c(gs0, gsSpline0),
                                                g1 = c(gs1, gsSpline1))
        } else {
            ## Now generate the vectors for Y * S(D, Z).

            if (!is.null(subset_index)) {
                newsubset <- subset_index
            } else {
                newsubset <- seq(1, nrow(data))
            }

            yvec <- as.vector(data[newsubset, yvar])
            dvec <- as.vector(data[newsubset, dvar])

            yvec <- yvec * (sest$sw1[, j] * dvec + sest$sw0[, j] * (1 - dvec))

            names(yvec) <- newsubset

            sset[[paste0("s", scount)]] <- list(beta = sest$beta[j],
                                                g0 = cbind(gs0, gsSpline0),
                                                g1 = cbind(gs1, gsSpline1),
                                                ys = yvec)
        }

        ## update counter (note scount is not referring
        ## to the list of IV regressions, but the components
        ## from the IV regressions)
        scount <- scount + 1
    }

    return(list(sset = sset, scount = scount))
}

#' GMM estimate of TE under point identification
#'
#' If the user sets the argument \code{point = TRUE} in the function
#' \code{ivmte}, then it is assumed that the treatment effect
#' parameter is point identified. The observational equivalence
#' condition is then set up as a GMM problem. Solving this GMM problem
#' recovers the coefficients on the MTR functions m0 and m1. Combining
#' these coefficients with the target gamma moments allows us to
#' estimate the target treatment effect.
#' @param sset a list of lists constructed from the function
#'     \link{gensset.mst}. Each inner list should include a
#'     coefficient corresponding to a term in an IV specification, a
#'     matrix of the estimates of the gamma moments conditional on (X,
#'     Z) for d = 0, and a matrix of the estimates of the gamma
#'     moments conditional on (X, Z) for d = 1. The column means of
#'     the last two matrices is what is used to generate the gamma
#'     moments.
#' @param gstar0 vector, the target gamma moments for d = 0.
#' @param gstar1 vector, the target gamma moments for d = 1.
#' @param itermax integer, maximum number of iterations allowed in the
#'     iterative GMM process. By default this is set to 2 (two-step
#'     GMM).
#' @param tol tolerance level for iterative GMM to terminate.
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @return a list containing the point estimate of the treatment
#'     effects, the standard errors, the 90% and 95% confidence
#'     intervals, the convergence code (see
#'     \code{\link[stats]{optim}}), the coefficients on the MTR, and
#'     the variance/covariance matrix of the MTR coefficient
#'     estimates.
gmmEstimate <- function(sset, gstar0, gstar1,
                        itermax = 2, tol = 1e-08, noisy = TRUE) {

    gmmMat <- NULL
    yMat   <- NULL

    for (s in 1:length(sset)) {

        ids <- as.integer(rownames(sset[[s]]$g0))

        gmmAdd <- cbind(ids, s,
                        sset[[s]]$g0,
                        sset[[s]]$g1)

        gmmMat <- rbind(gmmMat, gmmAdd)

        yAdd <- cbind(ids, s, sset[[s]]$ys)
        yMat <- rbind(yMat, yAdd)
    }

    N <- length(ids)

    gmmMat <- gmmMat[order(gmmMat[, 1], gmmMat[, 2]), ]
    yMat   <- yMat[order(yMat[, 1], yMat[, 2]), ]

    ids <- unique(gmmMat[, 1])
    gmmCompN <- ncol(gmmMat) - 2

    if (gmmCompN > length(sset)) {
        stop(gsub("\\s+", " ",
                  paste0("System is underidentified: excluding
                         target moments, there are ",
                         gmmCompN,
                         " unknown parameters/MTR coefficients and ",
                         length(sset),
                         " moment conditions (defined by IV-like
                         specifications). Either expand the number of
                         IV-like specifications, or modify m0 and m1.")))
    }

    gmmMat <- gmmMat[, -c(1, 2)]
    yMat   <- yMat[, -c(1, 2)]

    ## Perform iterative estimation
    theta <- rep(0, ncol(gmmMat))
    i <- 1
    diff <- Inf

    if (itermax > 2) warning("Itermax is capped at 2.")

    ## itermax is capped at 2, although it can be increased to
    ## correspond to iterated FGLS
    while (i <= itermax & i <= 2 & diff > tol) {

        if (i == 1) {
            thetaNew <- solve(t(gmmMat) %*% gmmMat) %*% t(gmmMat) %*% yMat
        } else {

            olsA <- lapply(ids, function(x) {
                gmmi <- gmmMat[as.integer(rownames(gmmMat)) == x, ]
                t(gmmi) %*% ematInv %*% gmmi
            })
            olsA <- Reduce("+", olsA)

            olsB <- lapply(ids, function(x) {
                gmmi <- gmmMat[as.integer(rownames(gmmMat)) == x, ]
                yvec <- yMat[as.integer(rownames(errors)) == x]
                t(gmmi) %*% ematInv %*% yvec
            })
            olsB <- Reduce("+", olsB)

            thetaNew <- solve(olsA) %*% olsB
        }

        errors <- yMat - gmmMat %*% thetaNew

        if (i <= (itermax - 1)) {
            emat <- lapply(ids, function(x) {
                evec <- errors[as.integer(rownames(errors)) == x]
                ## evec <- round(evec, 8)
                evec %*% t(evec)
            })
            emat <- Reduce("+", emat) / N
            emat <- (emat + t(emat)) / 2

            ematInv <- solve(emat)
            ematInv <- (ematInv + t(ematInv)) / 2
        }

        diff <- sqrt(sum((thetaNew - theta) ^ 2))
        theta <- thetaNew

        i <- i + 1
    }

    ## Construct point estimate and CI of TE

    rownames(theta) <- c(paste0("m0.", colnames(gstar0)),
                         paste0("m1.", colnames(gstar1)))

    te <- sum(c(colMeans(gstar0), colMeans(gstar1)) * theta)

    if (noisy == TRUE) {
        message()
        message(paste0("Treatment effect: ", round(te, 4)))
    }
    return(list(te = as.numeric(te),
                coef = theta))
}
