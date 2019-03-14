utils::globalVariables("u")

#' Estimation procedure from Mogstad, Torgovitsky (2017)
#'
#' This function estimates bounds on treatment effect parameters,
#' following the procedure described in Mogstad, Torgovitsky
#' (2017). Of the  parameters, the user can choose from the ATE,
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
#' Alternatively, point estimates can be obtained. Standard errors for
#' point estimates can be constructed using the bootstrap.
#'
#' @import methods stats utils
#'
#' @param bootstraps integer, default set to 0.
#' @param bootstraps.m integer, default set to size of data
#'     set. Determines the size of the subsample drawn from the
#'     original data set when performing inference via the
#'     bootstrap. This option applies only to the case of constructing
#'     confidence intervals for treatment effect bounds, i.e. it does
#'     not apply when \code{point = TRUE}.
#' @param bootstraps.replace boolean, default set to \code{TRUE}. This
#'     determines whether the resampling procedure used for inference
#'     will sample with replacement.
#' @param levels vector, real numbers between 0 and 1. Values
#'     correspond to the level of the confidence intervals constructed
#'     via bootstrap.
#' @param ci.type character, default set to 'both'. Set to 'forward'
#'     to construct the forward confidence interval for the treatment
#'     effect bound. Set to 'backward' to construct the backward
#'     confidence interval for the treatment effect bound. Set to
#'     'both' to construct both types of confidence intervals.
#' @param pvalue.tol numeric, default set to 1e-08. Tolerance level
#'     for determining p-value of treatment effect bound.
#' @param ivlike formula or vector of formulas used to specify the
#'     regressions for the IV-like estimands.
#' @param data \code{data.frame} used to estimate the treatment
#'     effects.
#' @param subset single subset condition or list of subset conditions
#'     corresponding to each IV-like estimand. The input must be
#'     logical. See \code{\link{l}} on how to input the argument. If
#'     the user wishes to select specific rows, construct a binary
#'     variable in the data set, and set the condition to use only
#'     those observations for which the binary variable is 1, e.g. the
#'     binary variable is \code{use}, and the subset condition is
#'     \code{use == 1}.
#' @param components a list of vectors of the terms/components from
#'     the regressions specifications we want to include in the set of
#'     IV-like estimands. To select the intercept term, include in the
#'     vector of variable names, `intercept'. If the the factorized
#'     counterpart of a variable \code{x = 1, 2, 3} is included in the
#'     IV-like specifications via \code{factor(x)}, the user can
#'     select the coefficients for specific factors by declaring the
#'     components \code{factor(x)-1, factor(x)-2, factor(x)-3}. See
#'     \code{\link{l}} on how to input the argument. If no components
#'     for a IV specification are given, then all components from that
#'     IV specification will be included.
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
#'     equivalence other than statistical noise, and the assumption
#'     that the model is correctly specified.
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
#'     point identified. If set to \code{TRUE}, then a FGLS procedure
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
#'     iterations allowed for FGLS estimation under point
#'     identification. So default estimate is the two-step FGLS.
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
#' jvec <- l(d, d, d, d)
#' svec <- l(, , , z %in% c(2, 4))
#'
#' ivmte(ivlike = ivlikespecs,
#'       data = dtm,
#'       components = jvec,
#'       propensity = d ~ z,
#'       subset = svec,
#'       m0 = ~  u + I(u ^ 2),
#'       m1 = ~  u + I(u ^ 2),
#'       uname = u,
#'       target = "att",
#'       m0.dec = TRUE,
#'       m1.dec = TRUE,
#'       bootstraps = 0,
#'       lpsolver = "lpSolveAPI")
#'
#' @export
ivmte <- function(bootstraps = 0, bootstraps.m,
                  bootstraps.replace = TRUE,
                  levels = c(0.99, 0.95, 0.90), ci.type = 'both',
                  pvalue.tol = 1e-08, ivlike, data, subset,
                  components, propensity, link, treat, m0, m1,
                  uname = u, target, target.weight0 = NULL,
                  target.weight1, target.knots0,
                  target.knots1 = NULL, late.Z, late.from, late.to,
                  late.X, eval.X, genlate.lb, genlate.ub,
                  obseq.tol = 0.05, grid.nu = 20, grid.nx = 20,
                  audit.nx = 2, audit.nu = 3, audit.max = 5,
                  audit.tol = 1e-08, m1.ub, m0.ub, m1.lb, m0.lb,
                  mte.ub, mte.lb, m0.dec, m0.inc, m1.dec, m1.inc,
                  mte.dec, mte.inc, lpsolver = NULL, point = FALSE,
                  point.itermax = 2, point.tol = 1e-08, noisy = TRUE) {

    call <- match.call(expand.dots = FALSE)

    ##---------------------------
    ## 1. Check linear programming dependencies
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
    ## 2. Check format of `formula', `subset', and `component' inputs
    ##---------------------------

    if (classFormula(ivlike)) {
        ivlike <- c(ivlike)
    }
    
    if (classList(ivlike)) {

        ## Convert formula, components, and subset inputs into lists
        length_formula <- length(ivlike)

        userComponents <- FALSE
        if (hasArg(components)) {
            if (!is.null(components)) {
                userComponents <- TRUE
            }
        }

        if (userComponents) {
            if (classList(components)) {
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
                components <- compList
            }
        } else {
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
            subset <- list(substitute(subset))
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
    ## 3. Check numeric arguments and case completion
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

        target.weight0 <- NULL
        target.weight1 <- NULL

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
        stop("Cannot set 'obseq.tol' below 0.")
    }
    if (!((grid.nu %% 1 == 0) & grid.nu >= 2)) {
        stop("grid.nu must be an integer greater than or equal to 2.")
    }
    if (!((grid.nx %% 1 == 0) & grid.nx >= 0)) {
        stop("'grid.nx' must be an integer greater than or equal to 0.")
    }

    if (!((audit.nx %% 1 == 0) & audit.nx > 0)) {
        stop("'audit.nx' must be an integer greater than or equal to 1.")
    }


    if (hasArg(m0.dec) | hasArg(m0.inc) |
        hasArg(m1.dec) | hasArg(m1.inc) |
        hasArg(mte.dec) | hasArg(mte.inc) |
        hasArg(m0.lb) | hasArg(m0.ub) |
        hasArg(m1.lb) | hasArg(m1.ub) |
        hasArg(mte.lb) | hasArg(mte.ub)) {

        noshape = FALSE ## indicator for whether shape restrictions declared

        if (!((audit.nu %% 1 == 0) & audit.nu > 0) | audit.nu < 2) {
            stop("'audit.nu' must be an integer greater than or equal to 2.")
        }

    } else {

        noshape = TRUE

        if (!((audit.nu %% 1 == 0) & audit.nu > 0)) {
            stop("'audit.nu' must be an integer greater than or equal to 1.")
        }
    }

    if (!((audit.max %% 1 == 0) & audit.max > 0)) {
        stop("'audit.max' must be an integer greater than or equal to 1.")
    }

    ## Bootstrap checks
    if (bootstraps < 0 | bootstraps %% 1 != 0) {
        stop("'bootstraps' must be an integer greater than or equal to 0.")
    }

    if (hasArg(bootstraps.m)) {
        if (bootstraps.m < 0 | bootstraps.m %% 1 != 0) {
            stop(gsub("\\s+", " ",
                      "'bootstraps.m' must be an integer greater than or equal
                      to 1."))
        }
    } else {
        bootstraps.m <- nrow(data)
    }

    if (!is.logical(bootstraps.replace))
        stop("'bootstraps.replace' must be TRUE or FALSE.")


    if (max(levels) >= 1 | min(levels) <= 0) {
        stop(gsub("\\s+", " ",
                  "'levels' must be a vector of values strictly between 0 and
                  1."))
    }
    levels <- sort(levels)


    if (! ci.type %in% c("forward", "backward", "both")) {
        stop(gsub("\\s+", " ",
                  "'ci.types' selects the type of confidence intervals to be
                  constructed for the treatment effect bound. It must be set to
                  either 'forward', 'backward', or 'both'."))
    }

    ##---------------------------
    ## 4. Restrict data to complete observations
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

    if (classFormula(ivlike)) {
        vars_formulas_x <- getXZ(ivlike)
        vars_formulas_z <- getXZ(ivlike, inst = TRUE)
        vars_y <- all.vars(ivlike)[1]
        ## terms_formulas <- attr(terms(Formula::as.Formula(ivlike)),
        ##                        "term.labels")
        terms_formulas_x <- getXZ(ivlike, terms = TRUE, inst = FALSE)
        terms_formulas_z <- getXZ(ivlike, terms = TRUE, inst = TRUE)
        
        
    } else if (classList(ivlike)) {
        if(!min(unlist(lapply(ivlike, classFormula)))) {
            stop(gsub("\\s+", " ",
                      "Not all elements in list of formulas are specified
                      correctly."))
        } else {
            vars_formulas_x <- unlist(lapply(ivlike,
                                             getXZ,
                                             inst = FALSE))
            vars_formulas_z <- unlist(lapply(ivlike,
                                             getXZ,
                                             inst = TRUE))

            vars_y <- unique(unlist(lapply(ivlike,
                                           function(x) all.vars(x)[[1]])))

            terms_formulas_x <- lapply(ivlike,
                                       getXZ,
                                       inst = FALSE,
                                       terms = TRUE)

            terms_formulas_z <- lapply(ivlike,
                                       getXZ,
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

        if (classList(subset)) {
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
    if (classFormula(m1) & classFormula(m0)) {
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
        if (classFormula(propensity)) {
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


            if (length(Formula::as.Formula(propensity))[1] == 0 &&
                length(all.vars(propensity)) > 1) {
                stop(gsub("\\s+", " ",
                          paste0("'propensity' argument must either be a
                          two-sided formula (if the propensity score is to be
                          estimated from the data), or a one-sided formula
                          containing a single variable on the RHS (where the
                          variable listed is included in the data, and
                          corresponds to propensity scores.")))
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

            propensity <- deparse(substitute(propensity))
            propensity <- Formula::as.Formula(paste("~", propensity))

            if (length(all.vars(propensity)) > 1) {
                stop(gsub("\\s+", " ",
                          paste0("'propensity' argument must either be a
                          two-sided formula (if the propensity score is to be
                          estimated from the data), or a one-sided formula
                          containing a single variable on the RHS (where the
                          variable listed is included in the data, and
                          corresponds to propensity scores.")))
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
                            paste(unique(c("intercept", terms_propensity)),
                                  collapse = " + "),
                            sep = " ~ ")

        propensity <- gsub("intercept", "1", propensity)

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

    if (classFormula(ivlike)) {
        comp_filler <- unstring(terms_formulas_x)
    } else {
        comp_filler <- lapply(terms_formulas_x,
                              function(x) as.character(unstring(x)))
    }

    ## Fill in components list if necessary
    if (userComponents) {
        compMissing1 <- unlist(lapply(components, function(x) {
            Reduce(paste, deparse(x)) == ""
        }))
        compMissing2 <- unlist(lapply(components, function(x) x == ""))
        compMissing <- as.logical(compMissing1 + compMissing2)

        if (sum(compMissing) > 0 & specCompWarn) {
            warning(gsub("\\s+", " ",
                         "Specifications without corresponding
                         component vectors will include all covariates when
                         constructing the S-set."),
                    call. = FALSE)
        }

        if (sum(compMissing) > 0) {
            components[compMissing] <- comp_filler[compMissing]
        }
    } else {
        print("am i doing the comp filler?")
        components <- comp_filler
    }

    print("componnets")
    print(components)
    ## stop("end of test")
    
    ## Keep only complete cases
    varError <- allvars[! allvars %in% colnames(data)]
    varError <- varError[varError != "intercept"]
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

    ##---------------------------
    ## 5. Implement estimates
    ##---------------------------

    estimateCall <- modcall(call,
                            newcall = ivmteEstimate,
                            dropargs = c("m0", "m1",
                                         "bootstraps", "data",
                                         "bootstraps.m",
                                         "bootstraps.replace",
                                         "subset",
                                         "levels", "ci.type",
                                         "treat", "propensity",
                                         "components", "lpsolver",
                                         "target.weight0", "target.weight1",
                                         "target.knots0", "target.knots1"),
                            newargs = list(m0 = quote(m0),
                                           m1 = quote(m1),
                                           target.weight0 =
                                               quote(target.weight0),
                                           target.weight1 =
                                               quote(target.weight1),
                                           target.knots0 = quote(target.knots0),
                                           target.knots1 = quote(target.knots1),
                                           data = quote(data),
                                           subset = quote(subset),
                                           lpsolver = quote(lpsolver),
                                           vars_y = quote(vars_y),
                                           vars_mtr = quote(vars_mtr),
                                           terms_mtr0 = quote(terms_mtr0),
                                           terms_mtr1 = quote(terms_mtr1),
                                           treat = quote(treat),
                                           propensity = quote(propensity),
                                           splinesobj = quote(splinesobj),
                                           components = quote(components)))


    ## Estimate bounds
    if (point == FALSE) {

        origEstimate <- eval(estimateCall)

        if (abs(origEstimate$bound[2] - origEstimate$bound[1]) < point.tol &
            noshape == TRUE) {

            message(gsub("\\s+", " ",
                         paste0("Width of bounds is narrower than 'point.tol',
                                which is set to ", point.tol, ". Estimation
                                under point identification will instead be
                                performed.")))

            estimateCall <-
                modcall(call,
                        newcall = ivmteEstimate,
                        dropargs = c("m0", "m1",
                                     "bootstraps", "data", "point",
                                     "bootstraps.m",
                                     "bootstraps.replace", "subset",
                                     "levels",
                                     "ci.type", "treat", "propensity",
                                     "components", "lpsolver",
                                     "target.weight0", "target.weight1",
                                     "target.knots0", "target.knots0"),
                        newargs = list(m0 = quote(m0),
                                       m1 = quote(m1),
                                       target.weight0 =
                                           quote(target.weight0),
                                       target.weight1 =
                                           quote(target.weight1),
                                       target.knots0 = quote(target.knots0),
                                       target.knots1 = quote(target.knots1),
                                       data = quote(data),
                                       subset = quote(subset),
                                       lpsolver = quote(lpsolver),
                                       point = TRUE,
                                       vars_y = quote(vars_y),
                                       vars_mtr = quote(vars_mtr),
                                       terms_mtr0 = quote(terms_mtr0),
                                       terms_mtr1 = quote(terms_mtr1),
                                       treat = quote(treat),
                                       propensity = quote(propensity),
                                       splinesobj = quote(splinesobj),
                                       components = quote(components)))

            return(eval(estimateCall))
        } else {
            ## Estimate bounds without resampling
            if (bootstraps == 0) {
                return(origEstimate)
            }
            ## Estimate bounds with resampling
            if (bootstraps > 0) {
                boundEstimates <- NULL

                b <- 1
                bootFailN <- 0
                bootFailNote <- ""
                bootFailIndex <- NULL

                while (b <= bootstraps) {
                    bootIDs  <- sample(seq(1, nrow(data)),
                                    size = nrow(data),
                                    replace = TRUE)
                    bdata <- data[bootIDs, ]

                    bootCall <-
                        modcall(call,
                                newcall = ivmteEstimate,
                                dropargs = c("m0", "m1",
                                             "bootstraps", "data",
                                             "noisy", "bootstraps.m",
                                             "bootstraps.replace",
                                             "subset", "levels", "ci.type",
                                             "treat",
                                             "propensity", "components",
                                             "lpsolver",
                                             "target.weight0", "target.weight1",
                                             "target.knots0", "target.knots1"),
                                newargs = list(m0 = quote(m0),
                                               m1 = quote(m1),
                                               target.weight0 =
                                                   quote(target.weight0),
                                               target.weight1 =
                                                   quote(target.weight1),
                                               target.knots0 =
                                                   quote(target.knots0),
                                               target.knots1 =
                                                   quote(target.knots1),
                                               data = quote(bdata),
                                               subset = quote(subset),
                                               lpsolver = quote(lpsolver),
                                               noisy = FALSE,
                                               vars_y = quote(vars_y),
                                               vars_mtr = quote(vars_mtr),
                                               terms_mtr0 = quote(terms_mtr0),
                                               terms_mtr1 = quote(terms_mtr1),
                                               treat = quote(treat),
                                               propensity = quote(propensity),
                                               splinesobj = quote(splinesobj),
                                               components = quote(components)))

                    bootEstimate <- try(eval(bootCall), silent = TRUE)
                    if (is.list(bootEstimate)) {
                        boundEstimates  <- rbind(boundEstimates,
                                                 bootEstimate$bound)

                        if (noisy == TRUE) {
                            message(paste0("Bootstrap iteration ", b,
                                           bootFailNote, "..."))
                        }
                        b <- b + 1
                        bootFailN <- 0
                        bootFailNote <- ""
                    } else {

                        if (noisy == TRUE) {
                            message(paste0("Bootstrap iteration ", b,
                                           bootFailNote,
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
                                        " failed. Failed bootstraps are
                                        repeated.")))
                }

                ## Obtain standard errors of bounds
                bootSE <- apply(boundEstimates, 2, sd)

                ## Construct confidence intervals
                if (ci.type == "backward" | ci.type == "forward") {
                    ci <- boundCI(bound = origEstimate$bound,
                                  bound.resamples = boundEstimates,
                                  n = nrow(data),
                                  m = bootstraps.m,
                                  levels = levels,
                                  type = ci.type)
                }

                if (ci.type == "both") {
                    ci <- list()
                    ci$backward <- boundCI(bound = origEstimate$bound,
                                           bound.resamples = boundEstimates,
                                           n = nrow(data),
                                           m = bootstraps.m,
                                           levels = levels,
                                           type = "backward")

                    ci$forward <- boundCI(bound = origEstimate$bound,
                                          bound.resamples = boundEstimates,
                                          n = nrow(data),
                                          m = bootstraps.m,
                                          levels = levels,
                                          type = "forward")
                }

                ci <- boundCI(bound = origEstimate$bound,
                              bound.resamples = boundEstimates,
                              n = nrow(data),
                              m = bootstraps.m,
                              levels = levels,
                              type = ci.type)

                ## Obtain p-value
                if (ci.type == "backward") {
                    pvalue <- boundPValue(ci = ci,
                                          bound = origEstimate$bound,
                                          bound.resamples = boundEstimates,
                                          n = nrow(data),
                                          m = bootstraps.m,
                                          levels = levels,
                                          type = "backward",
                                          tol = pvalue.tol)
                    names(pvalue) <- "backward"
                }

                if (ci.type == "forward") {
                    pvalue <- boundPValue(ci = ci,
                                          bound = origEstimate$bound,
                                          bound.resamples = boundEstimates,
                                          n = nrow(data),
                                          m = bootstraps.m,
                                          levels = levels,
                                          type = "forward",
                                          tol = pvalue.tol)
                    names(pvalue) <- "forward"
                }

                if (ci.type == "both") {
                    pvalue <- c(boundPValue(ci = ci$backward,
                                            bound = origEstimate$bound,
                                            bound.resamples = boundEstimates,
                                            n = nrow(data),
                                            m = bootstraps.m,
                                            levels = levels,
                                            type = "backward",
                                            tol = pvalue.tol),
                                boundPValue(ci = ci$forward,
                                            bound = origEstimate$bound,
                                            bound.resamples = boundEstimates,
                                            n = nrow(data),
                                            m = bootstraps.m,
                                            levels = levels,
                                            type = "forward",
                                            tol = pvalue.tol))
                    names(pvalue) <- c("backward", "forward")
                }

                ## Return output
                return(c(origEstimate,
                         list(bound.se = bootSE,
                              bound.bootstraps = boundEstimates,
                              ci = ci,
                              pvalue = pvalue,
                              bootstraps = bootstraps,
                              failed.bootstraps = length(bootFailIndex))))
            }
        }
    }

    ## Point estimate  without resampling
    if (point == TRUE & bootstraps == 0) {
        return(eval(estimateCall))
    }

    ## Point estimate with resampling
    if (point == TRUE & bootstraps > 0) {

        origEstimate <- eval(estimateCall)

        teEstimates  <- NULL
        mtrEstimates <- NULL
        propEstimates <- NULL

        b <- 1
        bootFailN <- 0
        bootFailNote <- ""
        bootFailIndex <- NULL

        while (b <= bootstraps) {
            bootIDs  <- sample(seq(1, nrow(data)),
                                 size = nrow(data),
                                 replace = TRUE)
            bdata <- data[bootIDs, ]

            bootCall <- modcall(call,
                                newcall = ivmteEstimate,
                                dropargs = c("m0", "m1",
                                             "bootstraps", "data",
                                             "noisy", "treat",
                                             "subset",
                                             "propensity",
                                             "components",
                                             "lpsolver",
                                             "target.weight0", "target.weight1",
                                             "target.knots0", "target.knots1"),
                                newargs = list(m0 = quote(m0),
                                               m1 = quote(m1),
                                               target.weight0 =
                                                   quote(target.weight0),
                                               target.weight1 =
                                                   quote(target.weight1),
                                               target.knots0 =
                                                   quote(target.knots0),
                                               target.knots1 =
                                                   quote(target.knots1),
                                               data = quote(bdata),
                                               subset = quote(subset),
                                               lpsolver = quote(lpsolver),
                                               noisy = FALSE,
                                               vars_y = quote(vars_y),
                                               vars_mtr = quote(vars_mtr),
                                               terms_mtr0 = quote(terms_mtr0),
                                               terms_mtr1 = quote(terms_mtr1),
                                               treat = quote(treat),
                                               propensity = quote(propensity),
                                               splinesobj = quote(splinesobj),
                                               components = quote(components)))

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

        ## Conf. int. 2: percentile method, using Z statistics
        ci290 <- origEstimate$te + c(qnorm(0.05), qnorm(0.95)) * bootSE
        ci295 <- origEstimate$te + c(qnorm(0.025), qnorm(0.975)) * bootSE
        names(ci290) <- c("5%", "95%")
        names(ci295) <- c("2.5%", "97.5%")

        mtrci290 <- sweep(x = tcrossprod(c(qnorm(0.05), qnorm(0.95)), mtrSE),
                          MARGIN = 2, origEstimate$mtr.coef, FUN = "+")
        mtrci295 <- sweep(x = tcrossprod(c(qnorm(0.025), qnorm(0.975)), mtrSE),
                          MARGIN = 2, origEstimate$mtr.coef, FUN = "+")

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
                      prop.ci2.95 = t(propci295),
                      bootstraps = bootstraps,
                      failed.bootstraps = length(bootFailIndex))))
    }
}

#' Construct confidence intervals for treatment effects under partial
#' identification
#'
#' This function constructs the forward and backward confidence
#' intervals for the treatment effect under partial identification.
#'
#' @param bound vector, bound of the treatment effects under partial
#'     identification.
#' @param bound.resamples matrix, stacked bounds of the treatment
#'     effects under partial identification. Each row corresponds to a
#'     subset resampled from the original data set.
#' @param n integer, size of original data set.
#' @param m integer, size of resampled data sets.
#' @param levels vector, real numbers between 0 and 1. Values
#'     correspond to the level of the confidence intervals constructed
#'     via bootstrap.
#' @param type character. Set to 'forward' to construct the forward
#'     confidence interval for the treatment effect bound. Set to
#'     'backward' to construct the backward confidence interval for
#'     the treatment effect bound. Set to 'both' to construct both
#'     types of confidence intervals.
#' @return if \code{type} is 'forward' or 'backward', then the
#'     corresponding type of confidence interval for each level is
#'     returned. The output is in the form of a matrix, with each row
#'     corresponding to a level. If \code{type} is 'both', then a list
#'     is returned. One element of the list is the matrix of backward
#'     confidence intervals, and the other element of the list is the
#'     matrix of forward confidence intervals.
boundCI <- function(bound, bound.resamples, n, m, levels, type) {

    if (type == "both") output <- list()

    ## Rescale and center bounds

    boundLBmod <- sqrt(m) *
        (bound.resamples[, 1] - bound[1])
    boundUBmod <- sqrt(m) *
        (bound.resamples[, 2] - bound[2])

    ## Construct backward confidence interval
    if (type %in% c('backward', 'both')) {
        backwardCI <- cbind(bound[1] +
                            (1 / sqrt(n)) *
                            quantile(boundLBmod,
                                     0.5 * (1 - levels),
                                     type = 1),
                            bound[2] +
                            (1 / sqrt(n)) *
                            quantile(boundUBmod,
                                     0.5 * (1 + levels),
                                     type = 1))

        colnames(backwardCI) <- c("lb.backward", "ub.backward")
        rownames(backwardCI) <- levels

        if (type == "both") output$backward <- backwardCI
        if (type == "backward") output <- backwardCI
    }

    ## Construct forward confidence interval
    if (type %in% c('forward', 'both')) {
        forwardCI <- cbind(bound[1] -
                           (1 / sqrt(n)) *
                           quantile(boundLBmod,
                                    0.5 * (1 + levels),
                                    type = 1),
                           bound[2] -
                           (1 / sqrt(n)) *
                           quantile(boundUBmod,
                                    0.5 * (1 - levels),
                                    type = 1))

        colnames(forwardCI) <- c("lb.forward", "ub.forward")
        rownames(forwardCI) <- levels

        if (type == "both") output$forward <- forwardCI
        if (type == "forward") output <- forwardCI
    }

    return(output)
}

#' Construct p-values for treatment effects under partial
#' identification
#'
#' This function estimates the p-value for the treatment effect under
#' partial identification. p-values corresponding to forward and
#' backward confidence intervals can be returned.
#'
#' @param ci matrix or list. If \code{type} is set to 'forward' or
#'     'backward', then \code{ci} should be a matrix of forward or
#'     backward confidence intervals corresponding to the levels
#'     declared in the option \code{levels}. If \code{type} is set to
#'     'both', then \code{ci} should be a list of two elements. One
#'     element is a matrix of forward confidence intervals, and the
#'     other element is a matrix of backward confidence intervals.
#' @param bound vector, bound of the treatment effects under partial
#'     identification.
#' @param bound.resamples matrix, stacked bounds of the treatment
#'     effects under partial identification. Each row corresponds to a
#'     subset resampled from the original data set.
#' @param n integer, size of original data set.
#' @param m integer, size of resampled data sets.
#' @param levels vector, real numbers between 0 and 1. Values
#'     correspond to the level of the confidence intervals constructed
#'     via bootstrap.
#' @param type character. Set to 'forward' to construct the forward
#'     confidence interval for the treatment effect bound. Set to
#'     'backward' to construct the backward confidence interval for
#'     the treatment effect bound. Set to 'both' to construct both
#'     types of confidence intervals.
#' @param tol numeric, default set to 1e-08. The p-value is
#'     constructed by iteratively adjusting the confidence level to
#'     find a confidence interval that does not contain 0. When the
#'     adjustment of the confidence level falls below \code{tol}, no
#'     further iterations are performed.
#' @return If \code{type} is 'forward' or 'backward', a scalar p-value
#'     corresponding to the type of confidence interval is
#'     returned. If \code{type} is 'both', a vector of p-values
#'     corresponding to the forward and backward confidence intervals
#'     is returned.
boundPValue <- function(ci, bound, bound.resamples, n, m, levels,
                        type, tol = 1e-08) {

    ci <- rbind(c(bound[1], bound[1]), ci, c(-Inf, Inf))
    rownames(ci) <- c(0, levels, 1)

    inVec <- apply(ci, 1, function(x) 0 > x[1] & 0 < x[2])

    lbPos <- which(sapply(seq(1, length(inVec) - 1), function(x) {
        inVec[x] == inVec[x + 1]
    }) == FALSE)

    levelLB <- c(0, levels, 1)[lbPos]
    levelUB <- c(0, levels, 1)[lbPos + 1]

    while (levelUB - levelLB > tol) {

        midpoint <- 0.5 * (levelLB + levelUB)

        newCI <- boundCI(bound, bound.resamples, n = 1000, m = 1000,
                         levels = midpoint, type = "backward")

        if (0 > newCI[1] & 0 < newCI[2]) {
            levelUB <- midpoint
        } else {
            levelLB <- midpoint
        }
    }
    return(1 - levelLB)
}

#' Single iteration of estimation procedure from Mogstad, Torgovitsky,
#' Santos (2018)
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
#'     corresponding to each IV-like estimand. The input must be
#'     logical. See \code{\link{l}} on how to input the argument. If
#'     the user wishes to select specific rows, construct a binary
#'     variable in the data set, and set the condition to use only
#'     those observations for which the binary variable is 1, e.g. the
#'     binary variable is \code{use}, and the subset condition is
#'     \code{use == 1}.
#' @param components a list of vectors of the terms/components from
#'     the regressions specifications we want to include in the set of
#'     IV-like estimands. To select the intercept term, include in the
#'     vector of variable names, `intercept'. See \code{\link{l}} on
#'     how to input the argument. If no components for a IV
#'     specification are given, then all components from that IV
#'     specification will be included.
#' @param propensity formula or variable name corresponding to
#'     propensity to take up treatment. If a formula is declared, then
#'     the function estimates propensity score according to the
#'     formula and link specified. If a variable name is declared,
#'     then the corresponding column in the data is taken as the
#'     vector of propensity scores.
#' @param link name of link function to estimate propensity score. Can
#'     be chosen from \code{linear}, \code{probit}, or \code{logit}.
#' @param treat variable name for treatment indicator.
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
#' @param vars_y character, variable name of observed outcome variable.
#' @param vars_mtr character, vector of variables entering into
#'     \code{m0} and \code{m1}.
#' @param terms_mtr0 character, vector of terms entering into
#'     \code{m0}.
#' @param terms_mtr1 character, vector of terms entering into
#'     \code{m1}.
#' @param splinesobj list of spline components in the MTRs for treated
#'     and control groups. Spline terms are extracted using
#'     \code{\link{removeSplines}}.
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
#'     equivalence other than statistical noise, and the assumption
#'     that the model is correctly specified.
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
#'     point identified. If set to \code{TRUE}, then a FGLS procedure
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
#'     iterations allowed for FGLS estimation under point
#'     identification. So default estimate is the two-step FGLS.
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @return Returns a list of results from throughout the estimation
#'     procedure. This includes all IV-like estimands; the propensity
#'     score model; bounds on the treatment effect; the estimated
#'     expectations of each term in the MTRs; the components and
#'     results of the LP problem.
ivmteEstimate <- function(ivlike, data, subset, components,
                           propensity, link, treat, m0, m1,
                           vars_y, vars_mtr, terms_mtr0, terms_mtr1,
                           splinesobj, uname = u,
                           target, target.weight0, target.weight1,
                           target.knots0 = NULL, target.knots1 = NULL,
                           late.Z, late.from, late.to, late.X, eval.X,
                           genlate.lb, genlate.ub, obseq.tol = 0.05,
                           grid.nu = 20, grid.nx = 20, audit.nx = 2,
                           audit.nu = 3, audit.max = 5,
                           audit.tol = 1e-08, m1.ub, m0.ub, m1.lb,
                           m0.lb, mte.ub, mte.lb, m0.dec, m0.inc,
                           m1.dec, m1.inc, mte.dec, mte.inc,
                           lpsolver = NULL, point = FALSE,
                           point.itermax = 2, point.tol = 1e-08,
                           noisy = TRUE) {

    call <- match.call(expand.dots = FALSE)

    ##---------------------------
    ## 1. Obtain propensity scores
    ##---------------------------

    if (noisy == TRUE) {
        message("Obtaining propensity scores...\n")
    }

    ## Estimate propensity scores

    pcall <- modcall(call,
                     newcall = propensity,
                     keepargs = c("link", "late.Z", "late.X"),
                     dropargs = "propensity",
                     newargs = list(data = quote(data),
                                    formula = propensity))
    pmodel <- eval(pcall)

    ##---------------------------
    ## 2. Generate target moments/gamma terms
    ##---------------------------

    if (noisy == TRUE) {
        message("Generating target moments...\n")
    }

    ## Parse polynomials
    if (!is.null(m0)) {
        m0call <- modcall(call,
                          newcall = polyparse,
                          keepargs = c("uname"),
                          newargs = list(formula = m0,
                                         data = quote(data)))
        pm0 <- eval(as.call(m0call))
    } else {
        pm0 <- NULL
    }

    if (!is.null(m1)) {
        m1call <- modcall(call,
                          newcall = polyparse,
                          keepargs = c("uname"),
                          newargs = list(formula = m1,
                                         data = quote(data)))
        pm1 <- eval(as.call(m1call))
    } else {
        pm1 <- NULL
    }

    ## Generate target weights
    if (!hasArg(target.weight0) & !hasArg(target.weight1)) {
        target.weight0 <- NULL
        target.weight1 <- NULL
    }

    if (is.null(target.weight0) & is.null(target.weight1)) {
        gentargetcall <- modcall(call,
                                 newcall = genTarget,
                                 keepargs = c("treat", "m1", "m0",
                                              "target", "late.Z",
                                              "late.from", "late.to",
                                              "late.X", "eval.X",
                                              "genlate.lb",
                                              "genlate.ub",
                                              "uname", "splinesobj",
                                              "point", "noisy"),
                                 dropargs = "data",
                                 newargs = list(data = quote(data),
                                                pmodobj = pmodel,
                                                pm0 = quote(pm0),
                                                pm1 = quote(pm1)))

    } else {
        gentargetcall <- modcall(call,
                                 newcall = genTarget,
                                 keepargs = c("treat", "m1", "m0",
                                              "target.weight0",
                                              "target.weight1",
                                              "target.knots0", "target.knots1",
                                              "uname", "splinesobj",
                                              "point", "noisy"),
                                 dropargs = "data",
                                 newargs = list(data = quote(data),
                                                pmodobj = pmodel,
                                                pm0 = quote(pm0),
                                                pm1 = quote(pm1)))
    }

    targetGammas <- eval(gentargetcall)
    gstar0 <- targetGammas$gstar0
    gstar1 <- targetGammas$gstar1

    ##---------------------------
    ## 3. Generate moments/gamma terms for IV-like estimands
    ##---------------------------

    if (noisy == TRUE) {
        message("Generating IV-like moments...")
    }

    sset  <- list() ## Contains all IV-like estimates and their
                    ## corresponding moments/gammas
    scount <- 1     ## counter for S-set constraints

    ## Construct `sset' object when a single IV-like specification is
    ## provided
    if (classFormula(ivlike)) {

        if (hasArg(subset)) {
            subset <- eval(subset[[1]], data)

            scall <- modcall(call,
                             newcall = ivEstimate,
                             keepargs = c("data", "components", "treat"),
                             newargs = list(formula = ivlike,
                                            subset = subset))

            sest <- eval(scall)
        } else {
            subset <- replicate(nrow(data), TRUE)
            scall <- modcall(call,
                             newcall = ivEstimate,
                             keepargs = c("data", "components", "treat"),
                             newargs = list(formula = ivlike,
                                            subset = subset))

            sest <- eval(scall)
        }


        ncomponents <- length(sest$betas)

        if (hasArg(subset)) {
            subset_index <- rownames(data[eval(substitute(subset), data), ])
        } else {
            subset_index <- rownames(data)
        }

        ## Generate moments (gammas) corresponding to IV-like
        ## estimands
        if (point == FALSE) {
            setobj <- genSSet(data = data,
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
            setobj <- genSSet(data = data,
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

    } else if (classList(ivlike)) {
        ## Construct `sset' object when multiple IV-like
        ## specifications are provided

        ## loop across IV specifications
        ivlikeCounter <- 1
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
            sest  <- ivEstimate(formula = sformula,
                                data = sdata,
                                components = scomponent,
                                treat = treat,
                                list = TRUE,
                                order = ivlikeCounter)
            ivlikeCounter <- ivlikeCounter + 1

            ## Generate moments (gammas) corresponding to IV-like
            ## estimands
            subset_index <- rownames(sdata)
            ncomponents <- length(sest$betas)
            pmodobj <- pmodel$phat[subset_index]
            if (point == FALSE) {
                setobj <- genSSet(data = data,
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
                setobj <- genSSet(data = data,
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

    ## Prepare FGLS estimate estimate if `point' agument is set to TRUE
    if (point == TRUE) {

        ## Obtain FGLS estimate

        fglsResult <- fglsEstimate(sset = sset,
                                 gstar0 = gstar0,
                                 gstar1 = gstar1,
                                 itermax = point.itermax,
                                 tol = point.tol,
                                 noisy = noisy)

        return(list(sset  = sset,
                    gstar = list(g0 = gstar0,
                                 g1 = gstar1),
                    propensity = pmodel,
                    te = fglsResult$te,
                    mtr.coef = fglsResult$coef))
    }

    ##---------------------------
    ## 4. Define constraint matrices using the audit
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
                          newcall = audit,
                          keepargs = audit.args,
                          newargs = list(data = quote(data),
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
        maxy <- max(data[, vars_y])
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
        miny <- min(data[, vars_y])
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
    ## 5. Obtain the bounds
    ##---------------------------

    audit <- eval(audit_call)

    message(paste0("Bound: (", audit$min, ", ", audit$max, ")\n"))

    ## include additional output material
    return(list(sset  = sset,
                gstar = list(g0 = gstar0,
                             g1 = gstar1),
                propensity = pmodel,
                bound = c(audit$min, audit$max),
                lpresult =  audit$lpresult,
                ## poly0 = pm0,
                ## poly1 = pm1,
                auditgrid = audit$gridobj,
                minobseq = audit$minobseq,
                splinesdict = list(splinesobj[[1]]$splinesdict,
                                   splinesobj[[2]]$splinesdict)))
}



#' Generating LP moments for IV-like estimands
#'
#' This function takes in the IV estimate and its IV-like
#' specification, and generates a list containing the corresponding
#' point estimate, and the corresponding moments (gammas) that will
#' enter into the constraint matrix of the LP problem.
#'
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
#' @param uname variable name for unobservable used in declaring MTRs
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
#' @param data \code{data.frame} used to estimate the treatment
#'     effects.
#' @param splinesobj list of spline components in the MTRs for treated
#'     and control groups. Spline terms are extracted using
#'     \code{\link{removeSplines}}.
#' @param pmodobj A vector of propensity scores.
#' @param pm0 A list of the monomials in the MTR for d = 0.
#' @param pm1 A list of the monomials in the MTR for d = 1.
#' @param point boolean, set to \code{FALSE} by default. \code{point}
#'     refers to whether the partial or point identification is
#'     desired. If set to \code{FALSE}, then the gamma moments are
#'     returned, i.e. sample averages are taken. If set to
#'     \code{TRUE}, then no sample averages are taken, and a matrix is
#'     returned. The sample average of each column of the matrix
#'     corresponds to a particular gamma moment.
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @return A list containing either the vectors of gamma moments for
#'     \code{D = 0} and \code{D = 1}, or a matrix of individual gamma
#'     values for \code{D = 0} and \code{D = 1}. Additoinally, two
#'     vectors are returned. \code{xindex0} and \code{xindex1} list
#'     the variables that interact with the unobservable \code{u} in
#'     \code{m0} and \code{m1}. \code{uexporder0} and
#'     \code{uexporder1} lists the exponents of the unobservable
#'     \code{u} in each term it appears in.
#'
#' @examples
#' ## Declare MTR functions
#' formula1 = ~ 1 + u
#' formula0 = ~ 1 + u
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Declare propensity score model
#' propensityObj <- propensity(formula = d ~ z,
#'                             data = dtm,
#'                             link = "linear")
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' polynomials1 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' ## Generate target gamma moments
#' genTarget(treat = "d",
#'           m0 = ~ 1 + u,
#'           m1 = ~ 1 + u,
#'           uname = u,
#'           target = "atu",
#'           data = dtm,
#'           splinesobj = splinesList,
#'           pmodobj = propensityObj,
#'           pm0 = polynomials0,
#'           pm1 = polynomials1,
#'           point = FALSE)
#'
#'
#' @export
genTarget <- function(treat, m0, m1, uname, target,
                      target.weight0, target.weight1,
                      target.knots0, target.knots1,
                      late.Z, late.from, late.to, late.X,
                      eval.X, genlate.lb, genlate.ub,
                      data, splinesobj, pmodobj, pm0, pm1,
                      point = FALSE, noisy = TRUE) {

    xindex0 <- NULL
    xindex1 <- NULL
    uexporder0 <- NULL
    uexporder1 <- NULL

    if (hasArg(target)) {
        if (target == "ate") {
            w1 <- wate1(data)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "att") {
            w1 <- watt1(data, mean(data[[treat]]), pmodobj$phat)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "atu") {
            w1 <- watu1(data, 1 - mean(data[[treat]]), pmodobj$phat)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "late") {
            if (!hasArg(late.X)) {
                late.X <- NULL
                eval.X <- NULL
            }
            w1 <- wlate1(data, late.from, late.to, substitute(late.Z),
                         pmodobj$model, substitute(late.X), eval.X)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "genlate") {
            w1 <- wgenlate1(data, genlate.lb, genlate.ub)
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
            if (point == FALSE) {
                gstar0 <- genGamma(pm0, w0$lb, w0$ub, w0$mp)
            } else {
                gstar0 <- genGamma(pm0, w0$lb, w0$ub, w0$mp, means = FALSE)
                xindex0 <- c(xindex0, pm0$xindex)
                uexporder0 <- c(uexporder0, pm0$exporder)
            }
        } else {
            gstar0 <- NULL
        }

        if (!is.null(m1)) {
            if (noisy == TRUE) {
                message("    Integrating terms for treated group...\n")
            }
            if (point == FALSE) {
                gstar1 <- genGamma(pm1, w1$lb, w1$ub, w1$mp)
            } else {
                gstar1 <- genGamma(pm1, w1$lb, w1$ub, w1$mp, means = FALSE)
                xindex1 <- c(xindex1, pm1$xindex)
                uexporder1 <- c(uexporder1, pm1$exporder)
            }
        } else {
            gstar1 <- NULL
        }

        if (point == FALSE) {
            gstarSplineObj0 <- genGammaSplines(splines = splinesobj[[1]],
                                               data = data,
                                               lb = w0$lb,
                                               ub = w0$ub,
                                               multiplier = w0$mp,
                                               d = 0)
            gstarSpline0 <- gstarSplineObj0$gamma

            gstarSplineObj1 <- genGammaSplines(splines = splinesobj[[2]],
                                               data = data,
                                               lb = w1$lb,
                                               ub = w1$ub,
                                               multiplier = w1$mp,
                                               d = 1)
            gstarSpline1 <- gstarSplineObj1$gamma
        } else {
            gstarSplineObj0 <- genGammaSplines(splines = splinesobj[[1]],
                                               data = data,
                                               lb = w0$lb,
                                               ub = w0$ub,
                                               multiplier = w0$mp,
                                               d = 0,
                                               means = FALSE)
            gstarSpline0 <- gstarSplineObj0$gamma
            xindex0 <- c(xindex0, gstarSplineObj0$interactions)
            uexporder0 <- c(uexporder0,
                            rep(-1, length(xindex0) - length(uexporder0)))

            gstarSplineObj1 <- genGammaSplines(splines = splinesobj[[2]],
                                               data = data,
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

                pm <- get(paste0("pm", d))

                gamma <- rep(0, length(pm$terms))

                for(i in 1:length(target.weight)) {

                    wKnotVarsL <- formalArgs(target.knots[[i]])
                    wKnotVarsU <- formalArgs(target.knots[[i + 1]])

                    if (wKnotVarsL[1] == "...") {
                        lb <- unlist(lapply(X = seq(1, nrow(data)),
                                            FUN = funEval,
                                            fun = target.knots[[i]]))
                    } else {
                        lb <- unlist(lapply(X = split(data[, wKnotVarsL],
                                                      seq(1, nrow(data))),
                                            FUN = funEval,
                                            fun = target.knots[[i]],
                                            argnames = wKnotVarsL))
                    }

                    if (wKnotVarsU[1] == "...") {
                        ub <- unlist(lapply(X = seq(1, nrow(data)),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]]))

                    } else {
                        ub <- unlist(lapply(X = split(data[, wKnotVarsU],
                                                      seq(1, nrow(data))),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]],
                                            argnames = wKnotVarsU))
                    }

                    wValVars  <- formalArgs(target.weight[[i]])

                    if (wValVars[1] == "...") {
                        weights <- unlist(lapply(X = seq(1, nrow(data)),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]]))
                    } else {
                        weights <- unlist(lapply(X = split(data[, wValVars],
                                                           seq(1, nrow(data))),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]],
                                                 argnames = wValVars))
                    }

                    if (point == FALSE) {
                        gamma <- gamma + genGamma(pm,
                                                  lb = lb,
                                                  ub = ub,
                                                  multiplier = weights)
                    } else {
                        gamma <- gamma + genGamma(pm,
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
                        lb <- unlist(lapply(X = seq(1, nrow(data)),
                                            FUN = funEval,
                                            fun = target.knots[[i]]))
                    } else {
                        lb <- unlist(lapply(X = split(data[, wKnotVarsL],
                                                      seq(1, nrow(data))),
                                            FUN = funEval,
                                            fun = target.knots[[i]],
                                            argnames = wKnotVarsL))
                    }

                    if (wKnotVarsU[1] == "...") {
                        ub <- unlist(lapply(X = seq(1, nrow(data)),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]]))

                    } else {
                        ub <- unlist(lapply(X = split(data[, wKnotVarsU],
                                                      seq(1, nrow(data))),
                                            FUN = funEval,
                                            fun = target.knots[[i + 1]],
                                            argnames = wKnotVarsU))
                    }

                    wValVars  <- formalArgs(target.weight[[i]])

                    if (wValVars[1] == "...") {
                        weights <- unlist(lapply(X = seq(1, nrow(data)),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]]))
                    } else {
                        weights <- unlist(lapply(X = split(data[, wValVars],
                                                           seq(1, nrow(data))),
                                                 FUN = funEval,
                                                 fun = target.weight[[i]],
                                                 argnames = wValVars))
                    }

                    if (point == FALSE) {

                        gammaSplines <- gammaSplines +
                            genGammaSplines(splines = noSplineMtr,
                                            data = data,
                                            lb = lb,
                                            ub = ub,
                                            multiplier = weights,
                                            d = d)$gamma
                    } else {
                        gammaSplines <- gammaSplines +
                            genGammaSplines(splines = noSplineMtr,
                                            data = data,
                                            lb = lb,
                                            ub = ub,
                                            multiplier = weights,
                                            d = d,
                                            means = FALSE)$gamma
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

    return(list(gstar0 = gstar0,
                gstar1 = gstar1,
                xindex0 = xindex0,
                xindex1 = xindex1,
                uexporder0 = uexporder0,
                uexporder1 = uexporder1))
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
#'
#' @examples
#' ## Declare empty list to be updated (in the event multiple IV like
#' ## specifications are provided)
#' sSet <- list()
#'
#' ## Declare MTR formulas
#' formula1 = ~ 1 + u
#' formula0 = ~ 1 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' polynomials1 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' ## Generate propensity score model
#' propensityObj <- propensity(formula = d ~ z,
#'                             data = dtm,
#'                             link = "linear")
#'
#' ## Generate IV estimates
#' ivEstimates <- ivEstimate(formula = ey ~ d | z,
#'                           data = dtm,
#'                           components = l(d),
#'                           treat = d,
#'                           list = FALSE)
#'
#' ## Construct S-set, which contains the coefficients and weights
#' ## coresponding to various IV-like estimands
#' genSSet(data = dtm,
#'         sset = sSet,
#'         sest = ivEstimates,
#'         splinesobj = splinesList,
#'         pmodobj = propensityObj$phat,
#'         pm0 = polynomials0,
#'         pm1 = polynomials1,
#'         ncomponents = 1,
#'         scount = 1)
#'
#' @export
genSSet <- function(data, sset, sest, splinesobj, pmodobj, pm0, pm1,
                        ncomponents, scount, subset_index, means = TRUE,
                        yvar, dvar, noisy = TRUE) {

    if (!hasArg(subset_index)) subset_index <- NULL

    for (j in 1:ncomponents) {
        if (noisy == TRUE) {
            message(paste0("    Moment ", scount, "..."))
        }

        if (!is.null(pm0)) {
            if (means == TRUE) {
                gs0 <- genGamma(monomials = pm0,
                                    lb = pmodobj,
                                    ub = 1,
                                    multiplier = sest$sw0[, j],
                                    subset = subset_index)
            } else {
                gs0 <- genGamma(monomials = pm0,
                                    lb = pmodobj,
                                    ub = 1,
                                    multiplier = sest$sw0[, j],
                                    subset = subset_index,
                                    means = FALSE)
            }
        } else {
            gs0 <- NULL
        }

        if (!is.null(pm1)) {
            if (means == TRUE) {
                gs1 <- genGamma(monomials = pm1,
                                    lb = 0,
                                    ub = pmodobj,
                                    multiplier = sest$sw1[, j],
                                    subset = subset_index)
            } else {
                gs1 <- genGamma(monomials = pm1,
                                    lb = 0,
                                    ub = pmodobj,
                                    multiplier = sest$sw1[, j],
                                    subset = subset_index,
                                    means = FALSE)
            }
        } else {
            gs1 <- NULL
        }

        if (means == TRUE) {
            gsSpline0 <- genGammaSplines(splines = splinesobj[[1]],
                                             data = data,
                                             lb = pmodobj,
                                             ub = 1,
                                             multiplier = sest$sw0[, j],
                                             subset = subset_index,
                                             d = 0)$gamma

            gsSpline1 <- genGammaSplines(splines = splinesobj[[2]],
                                             data = data,
                                             lb = 0,
                                             ub = pmodobj,
                                             multiplier = sest$sw1[, j],
                                             subset = subset_index,
                                             d = 1)$gamma
        } else {
            gsSpline0 <- genGammaSplines(splines = splinesobj[[1]],
                                             data = data,
                                             lb = pmodobj,
                                             ub = 1,
                                             multiplier = sest$sw0[, j],
                                             subset = subset_index,
                                             d = 0,
                                             means = FALSE)$gamma

            gsSpline1 <- genGammaSplines(splines = splinesobj[[2]],
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

#' FGLS estimate of TE under point identification
#'
#' If the user sets the argument \code{point = TRUE} in the function
#' \code{ivmte}, then it is assumed that the treatment effect
#' parameter is point identified. The observational equivalence
#' condition is then set up as a FGLS problem. Solving this FGLS problem
#' recovers the coefficients on the MTR functions m0 and m1. Combining
#' these coefficients with the target gamma moments allows us to
#' estimate the target treatment effect.
#' @param sset a list of lists constructed from the function
#'     \link{genSSet}. Each inner list should include a
#'     coefficient corresponding to a term in an IV specification, a
#'     matrix of the estimates of the gamma moments conditional on (X,
#'     Z) for d = 0, and a matrix of the estimates of the gamma
#'     moments conditional on (X, Z) for d = 1. The column means of
#'     the last two matrices is what is used to generate the gamma
#'     moments.
#' @param gstar0 vector, the target gamma moments for d = 0.
#' @param gstar1 vector, the target gamma moments for d = 1.
#' @param itermax integer, maximum number of iterations allowed in the
#'     iterative FGLS process. By default this is set to 2 (two-step
#'     FGLS).
#' @param tol tolerance level for iterative FGLS to terminate.
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
#'
#' @examples
#' ## Declare empty list to be updated (in the event multiple IV like
#' ## specifications are provided
#' sSet <- list()
#'
#' ## Declare MTR formulas
#' formula1 = ~ 0 + u
#' formula0 = ~ 0 + u
#'
#' ## Construct object that separates out non-spline components of MTR
#' ## formulas from the spline components. The MTR functions are
#' ## obtained from this object by the function 'genSSet'.
#' splinesList = list(removeSplines(formula0), removeSplines(formula1))
#'
#' ## Construct MTR polynomials
#' polynomials0 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#' polynomials1 <- polyparse(formula = formula0,
#'                  data = dtm,
#'                  uname = u,
#'                  as.function = FALSE)
#'
#' ## Generate propensity score model
#' propensityObj <- propensity(formula = d ~ z,
#'                             data = dtm,
#'                             link = "linear")
#'
#' ## Generate IV estimates
#' ivEstimates <- ivEstimate(formula = ey ~ d | z,
#'                           data = dtm,
#'                           components = l(intercept, d),
#'                           treat = d,
#'                           list = FALSE)
#'
#' ## Generate target gamma moments
#' targetGamma <- genTarget(treat = "d",
#'                          m0 = ~ 1 + u,
#'                          m1 = ~ 1 + u,
#'                          uname = u,
#'                          target = "atu",
#'                          data = dtm,
#'                          splinesobj = splinesList,
#'                          pmodobj = propensityObj,
#'                          pm0 = polynomials0,
#'                          pm1 = polynomials1,
#'                          point = TRUE)
#'
#' ## Construct S-set. which contains the coefficients and weights
#' ## corresponding to various IV-like estimands
#' sSet <- genSSet(data = dtm,
#'                 sset = sSet,
#'                 sest = ivEstimates,
#'                 splinesobj = splinesList,
#'                 pmodobj = propensityObj$phat,
#'                 pm0 = polynomials0,
#'                 pm1 = polynomials1,
#'                 ncomponents = 2,
#'                 scount = 1,
#'                 yvar = "ey",
#'                 dvar = "d",
#'                 means = FALSE)
#'
#' ## Obtain point estimates using FGLS
#' fglsEstimate(sset = sSet$sset,
#'              gstar0 = targetGamma$gstar0,
#'              gstar1 = targetGamma$gstar1)
#'
#' @export
fglsEstimate <- function(sset, gstar0, gstar1,
                        itermax = 2, tol = 1e-08, noisy = TRUE) {

    fglsMat <- NULL
    yMat   <- NULL

    for (s in 1:length(sset)) {

        ids <- as.integer(rownames(sset[[s]]$g0))

        fglsAdd <- cbind(ids, s,
                        sset[[s]]$g0,
                        sset[[s]]$g1)

        fglsMat <- rbind(fglsMat, fglsAdd)

        yAdd <- cbind(ids, s, sset[[s]]$ys)
        yMat <- rbind(yMat, yAdd)
    }

    N <- length(ids)

    fglsMat <- fglsMat[order(fglsMat[, 1], fglsMat[, 2]), ]
    yMat   <- yMat[order(yMat[, 1], yMat[, 2]), ]

    ids <- unique(fglsMat[, 1])
    fglsCompN <- ncol(fglsMat) - 2

    if (fglsCompN > length(sset)) {
        stop(gsub("\\s+", " ",
                  paste0("System is underidentified: excluding
                         target moments, there are ",
                         fglsCompN,
                         " unknown parameters/MTR coefficients and ",
                         length(sset),
                         " moment conditions (defined by IV-like
                         specifications). Either expand the number of
                         IV-like specifications, or modify m0 and m1.")))
    }

    fglsMat <- fglsMat[, -c(1, 2)]
    yMat   <- yMat[, -c(1, 2)]

    ## Perform iterative estimation
    theta <- rep(0, ncol(fglsMat))
    i <- 1
    diff <- Inf

    if (itermax > 2) warning("Itermax is capped at 2.")

    ## itermax is capped at 2, although it can be increased to
    ## correspond to iterated FGLS
    while (i <= itermax & i <= 2 & diff > tol) {

        if (i == 1) {
            thetaNew <- solve(t(fglsMat) %*% fglsMat) %*% t(fglsMat) %*% yMat
        } else {

            olsA <- lapply(ids, function(x) {
                fglsi <- fglsMat[as.integer(rownames(fglsMat)) == x, ]
                t(fglsi) %*% ematInv %*% fglsi
            })
            olsA <- Reduce("+", olsA)

            olsB <- lapply(ids, function(x) {
                fglsi <- fglsMat[as.integer(rownames(fglsMat)) == x, ]
                yvec <- yMat[as.integer(rownames(errors)) == x]
                t(fglsi) %*% ematInv %*% yvec
            })
            olsB <- Reduce("+", olsB)

            thetaNew <- solve(olsA) %*% olsB
        }

        errors <- yMat - fglsMat %*% thetaNew

        if (i <= (itermax - 1)) {
            emat <- lapply(ids, function(x) {
                evec <- errors[as.integer(rownames(errors)) == x]

                evec <- round(evec, 8)
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
