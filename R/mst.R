utils::globalVariables("u")

#' Instrumental Variables: Extrapolation by Marginal Treatment Effects
#'
#' This function provides a general framework for using the marginal treatment
#' effect (MTE) to extrapolate. The model is the same binary treatment
#' instrumental variable (IV) model considered by
#' \href{https://doi.org/10.2307/2951620}{Imbens and Angrist (1994)} and
#' \href{https://doi.org/10.1111/j.1468-0262.2005.00594.x}{Heckman and Vytlacil
#' (2005)}. The framework on which this function is based was developed by
#' \href{https://doi.org/10.3982/ECTA15463}{Mogstad, Santos and Torgovitsky
#' (2018)}. See also the recent survey paper on extrapolation in IV models by
#' \href{https://doi.org/10.1146/annurev-economics-101617-041813}{Mogstad and
#' Torgovitsky (2018)}.
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
#'     be chosen from \code{linear}, \code{probit}, or
#'     \code{logit}. Default is set to "logit".
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
#'     variable names in \code{data}. If the knots are constant across
#'     all observations, then the user can instead submit the vector
#'     of knots instead of a function.
#' @param target.knots1 user-defined set of functions defining the
#'     knots associated with splines weights for the treated
#'     group. The arguments of the function should be variable names
#'     in \code{data}. If the knots are constant across all
#'     observations, then the user can instead submit the vector of
#'     knots instead of a function.
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
#' @param initgrid.nu number of evenly spread points in the interval
#'     [0, 1] of the unobservable u used to form the initial grid for
#'     imposing shape restrictions on the MTRs.
#' @param initgrid.nx number of evenly spread points of the covariates
#'     to use to form the initial grid for imposing shape restrictions
#'     on the MTRs.
#' @param audit.nx number of points on the covariates space to audit
#'     in each iteration of the audit procedure.
#' @param audit.nu number of points in the interval [0, 1],
#'     corresponding to the normalized value of the unobservable term,
#'     to audit in each iteration of the audit procedure.
#' @param audit.add maximum number of points to add to the grids for
#'     imposing each kind of shape constraint. So if there are 5
#'     different kinds of shape constraints, there can be at most
#'     \code{audit.add * 5} additional points added to the grid.
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
#'     point identified. If set to \code{TRUE}, then a two-step GMM
#'     procedure is implemented to estimate the treatment
#'     effects. Shape constraints on the MTRs will be ignored under
#'     point identification.
#' @param point.eyeweight boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} if GMM point estimate should use the identity
#'     weighting matrix (i.e. one-step GMM).
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @param seed integer, the seed that determines the random grid in
#'     the audit procedure.
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
                  components, propensity, link = "logit", treat, m0,
                  m1, uname = u, target, target.weight0 = NULL,
                  target.weight1, target.knots0, target.knots1 = NULL,
                  late.Z, late.from, late.to, late.X, eval.X,
                  genlate.lb, genlate.ub, obseq.tol = 0.05,
                  initgrid.nu = 10, initgrid.nx = 20, audit.nx = 2500,
                  audit.nu = 25, audit.add = 100, audit.max = 25,
                  audit.tol = 1e-08, m1.ub, m0.ub, m1.lb, m0.lb,
                  mte.ub, mte.lb, m0.dec, m0.inc, m1.dec, m1.inc,
                  mte.dec, mte.inc, lpsolver = NULL, point = FALSE,
                  point.eyeweight = FALSE, noisy = TRUE, seed = 12345) {

    ## Include into Roxygen documentation once document is published..
    ## A detailed description of the module and its features
    ## can be found in
    ## \href{https://a-torgovitsky.github.io/shea-torgovitsky.pdf}{Shea and
    ## Torgovitsky (2019)}.

    call <- match.call(expand.dots = FALSE)

    ##---------------------------
    ## 1. Check linear programming dependencies
    ##---------------------------

    if (hasArg(lpsolver)) lpsolver <- tolower(lpsolver)

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
                              "rcplex",
                              "cplexapi",
                              "lpsolve",
                              "lpsolveapi")) {
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
    ## 2. Check format of non-numeric arguments
    ##---------------------------

    ## Character arguments will be converted to lowercase
    if (hasArg(target))   target   <- tolower(target)
    if (hasArg(link))     link     <- tolower(link)
    if (hasArg(ci.type))  ci.type  <- tolower(ci.type)

    ## Convert ivlike formula into a list (i.e. a one-element list
    ## with one formula), which is a more robust framework
    if (classFormula(ivlike)) {
        ivlike <- c(ivlike)
    }

    ## Convert formula, components, and subset inputs into lists
    length_formula <- length(ivlike)

    userComponents <- FALSE
    if (hasArg(components)) {
        if (!is.null(components)) {
            if (length(components) > 1) {
                userComponents <- TRUE
            }
            if (length(components) == 1) {
                if (Reduce(paste, deparse(components)) != "list()"){
                    userComponents <- TRUE
                }
            }
        }
    }

    if (userComponents) {
        length_components <- length(components)
        if (length_formula == 1) {
            ## When a single formula is provided, then the list of
            ## components should be treated as a single vector of
            ## components. The way in which the user declares the
            ## components can be problematic. The function must figure
            ## out if the components list is entered directly, or as a
            ## variable.
            componentsTmp <- gsub("\\s+", " ",
                                  Reduce(paste,
                                         deparse(substitute(components))))
            if (substr(componentsTmp, 1, 2) == "l(") {
                components <- deparse(substitute(components))
                components <- gsub("\\s+", " ", Reduce(paste, components))
                if (substr(componentsTmp, 1, 4) != "l(c(") {
                    internals <- substr(components, 3,
                                        nchar(components) - 1)
                    charList <- unique(unlist(strsplit(x = internals,
                                                       split = "")))
                    if (length(charList) == 2 && all(charList == c(",", " "))) {
                        internals <- ""
                        warning(gsub("\\s+", " ",
                                     "No list of components provided.
                                      All covariates in each
                                      IV-like specification will be included
                                      when constructing each S-set."),
                                call. = FALSE)
                    }
                    components <- paste0("l(c(", internals, "))")
                }
                components <- eval(parse(text = components))
                length_components <- 1
            } else {
                components <- unlist(lapply(components, deparse))
                components <- gsub("\\s+", " ", Reduce(paste, components))
                if (substr(components, 1, 2) == "c(") {
                    components <- paste0("l(", components, ")")
                } else {
                    components <- paste0("l(c(", components, "))")
                }
                components <- eval(parse(text = components))
            }
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
        ## If subsetting is not a list, convert it to a list
        if (!(classList(subset))) {
            ## Check if character, if so, then may need to
            ## deparse. Also check if logical, in which case less
            ## needs to be done.
            subsetChar <- suppressWarnings(
                try(is.character(subset), silent = TRUE))
            subsetLogic <- suppressWarnings(
                try(is.logical(subset), silent = TRUE))
            if (subsetChar == TRUE) {
                stop(gsub("\\s+", " ",
                          "Subset conditions should be logical
                           expressions involving variable names from the
                           data set and logical operators."))
            } else if (subsetLogic == TRUE) {
                ## Currently prohobit logical vectors. Instead
                ## restrict to expressions only. Below is the code to
                ## revert to the case where logical vectors are
                ## allowed.
                ##
                ## subsetList <- list()
                ## subsetList[[1]] <- subset
                ## subset <- subsetList
                stop(gsub("\\s+", " ",
                          "Subset conditions should be logical
                           expressions involving variable names from the
                           data set and logical operators. Make sure the name of
                           the variable is not the same as that of another
                           object in the workspace. If so, rename the object
                           in the workspace, or the variable in the data set."))
            } else {
                subsetStrSubs <- Reduce(paste, deparse(substitute(subset)))
                subsetStr <- suppressWarnings(
                    try(Reduce(paste, deparse(subset)), silent = TRUE))
                if (class(subsetStr) == "try-error") { ## i.e. logical
                    subset <- paste0("l(", subsetStrSubs, ")")
                } else { ## i.e. variable, for looping
                    subset <- paste0("l(", subsetStr, ")")
                }
            }
            if (subsetLogic != TRUE) {
                subsetLogicOp <- 0
                for (operator in c("==", "<", ">", "%in%")) {
                    subsetLogicOp <- subsetLogicOp + grepl(operator, subset)
                }
                if (subsetLogicOp == 0) {
                    stop(gsub("\\s+", " ",
                              "Subset conditions should be logical
                           expressions involving variable names from the
                           data set and logical operators."))
                }
                subset <- eval(parse(text = subset))
            }
        }
        ## Fill in any missing subset slots
        if (length(subset) > 1 && length(subset) != length_formula) {
            stop(gsub("\\s+", " ",
                      "Number of subset conditions not equal to number of
                       IV specifications. Either declare a single subset
                       condition to be applied to all IV specifications; or
                       declare a list of subset conditions, one for each IV
                       specificaiton. An empty element in the list of subset
                       conditions corresponds to using the full sample."))
        }
        if (length(subset) == 1 && length_formula > 1) {
            subset <- rep(subset, length_formula)
        }
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
            if ((hasArg(late.X) & !hasArg(eval.X)) |
               !hasArg(late.X) & hasArg(eval.X)) {
                stop(gsub("\\s+", " ",
                          "If the target parameter is 'late', then either both
                          late.X and eval.X are specified, or neither are
                          specified."))
            }

            ## Check the LATE arguments are declared properly.
            if (classList(late.to)) late.to <- unlist(late.to)
            if (classList(late.from)) late.from <- unlist(late.from)
            zIsVec <- substr(deparse(substitute(late.Z)), 1, 2) == "c("
            zIsList <- substring(deparse(substitute(late.Z)), 1, 2) == "l("
            if (zIsVec) {
                late.Z <- restring(substitute(late.Z), substitute = FALSE)
            } else if (zIsList) {
                late.Z <- restring(substitute(late.Z), substitute = FALSE,
                                   command = "l")
            } else {
                late.Z <- deparse(substitute(late.Z))
            }

            if (length(late.to) != length(late.from) |
                length(late.to) != length(late.Z)) {
                stop(gsub("\\s+", " ",
                          "The number of variables declared in 'late.Z'
                           must be equal to the length of the vectors
                           declared in 'late.to' and 'late.from'."))
            }
            if (!is.numeric(late.to) | !is.numeric(late.from)) {
                stop("Vectors 'late.from' and 'late.to' must be numeric.")
            }
            if (all(late.to == late.from)) {
                stop(gsub("\\s+", " ",
                          "'late.to' must be different from 'late.from'."))
            }

            ## Check that included insturments are declared properly
            if (hasArg(eval.X)) {
                if (classList(eval.X)) eval.X <- unlist(eval.X)
                xIsVec <- substr(deparse(substitute(late.X)), 1, 2) == "c("
                xIsList <- substring(deparse(substitute(late.X)), 1, 2) == "l("
                if (xIsVec) {
                    late.X <- restring(substitute(late.X),
                                         substitute = FALSE)
                } else if (xIsList) {
                    late.X <- restring(substitute(late.X),
                                         substitute = FALSE,
                                         command = "l")
                } else {
                    late.X <- deparse(substitute(late.X))
                }
                if (length(late.X) != length(eval.X)) {
                    stop(gsub("\\s+", " ",
                              "The number of variables declared in 'late.X'
                       must be equal to the length of the vector
                       declared in 'eval.X'."))
                }
                if (!is.numeric(eval.X)) {
                    stop("Vector 'eval.X' must be numeric.")
                }
            }
        }
        if (target == "genlate") {
            if (genlate.lb < 0 | genlate.ub > 1) {
                stop(gsub("\\s+", " ",
                          "'genlate.lb' and 'genlate.ub' must be between 0
                          and 1."))
            }
            if (genlate.lb >= genlate.ub) {
                stop("'genlate.lb' must be strictly less than 'genlate.ub'.")
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
                hasArg(audit.nu) | hasArg(audit.nx) | hasArg(audit.add) |
                hasArg(initgrid.nu) | hasArg(initgrid.nx)|
                hasArg(audit.tol) | hasArg(audit.max)) {
                warning(gsub("\\s+", " ",
                             "If argument 'point' is set to TRUE, then shape
                             restrictions on m0 and m1 are ignored, and the
                             audit procedure is not implemented."))
            }
        }
    }

    ## Seed check
    if (!(is.numeric(seed) & length(seed) == 1)) {
        stop("'seed' must be a scalar.")
    }

    ## Audit checks
    if (!(is.numeric(obseq.tol) & obseq.tol >= 0)) {
        stop("Cannot set 'obseq.tol' below 0.")
    }
    if (!((initgrid.nu %% 1 == 0) & initgrid.nu >= 2)) {
        stop("initgrid.nu must be an integer greater than or equal to 2.")
    }
    if (!((initgrid.nx %% 1 == 0) & initgrid.nx >= 0)) {
        stop("'initgrid.nx' must be an integer greater than or equal to 0.")
    }

    if (!((audit.nx %% 1 == 0) & audit.nx > 0)) {
        stop("'audit.nx' must be an integer greater than or equal to 1.")
    }
    if (!((audit.add %% 1 == 0) & audit.add > 0)) {
        stop("'audit.add' must be an integer greater than or equal to 1.")
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
    if (bootstraps < 0 | bootstraps %% 1 != 0 | bootstraps == 1) {
        stop(gsub("\\s+", " ",
                  "'bootstraps' must either be 0, or be an integer greater than
                   or equal to 2."))
    }

    if (hasArg(bootstraps.m)) {
        if (bootstraps.m < 0 | bootstraps.m %% 1 != 0) {
            stop(gsub("\\s+", " ",
                      "'bootstraps.m' must be an integer greater than or equal
                      to 1."))
        }
        if (point == TRUE) {
            warning(gsub("\\s+", " ",
                         "Argument 'bootstrap.m' is only used for
                          partial identification, and will be ignored
                          under point identification."), call. = FALSE)
        }
    } else {
        bootstraps.m <- nrow(data)
    }
    if (!is.logical(bootstraps.replace)) {
        stop("'bootstraps.replace' must be TRUE or FALSE.")
    }
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

    if (hasArg(point.eyeweight) && point == FALSE) {
        warning(gsub("\\s+", " ",
                     "Argument 'point.eyeweight' is only used for
                      point identification, and will be ignored when
                      point = FALSE."), call. = FALSE)
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
    vars_components <- c()

    terms_formulas_x <- c()
    terms_formulas_z <- c()
    terms_mtr0       <- c()
    terms_mtr1       <- c()
    terms_components <- c()

    if (classList(ivlike)) {

        if (!min(unlist(lapply(ivlike, classFormula)))) {
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
            svec <- gsub("\\s+", " ",
                         Reduce(paste, deparse(substitute(subset))))
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
        origm0 <- m0
        origm1 <- m1
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
            if (length(propensity) == 3) {
                ptreat <- all.vars(propensity)[1]
                vars_propensity <- all.vars(propensity)

                if (hasArg(treat)) {
                    if (ptreat != deparse(substitute(treat))) {
                        stop(gsub("\\s+", " ",
                                  "Variable listed in 'treat' argument
                                 differs from dependent variable in propensity
                                 score formula. Dependent variable from
                                 propensity score formula will be used as the
                                 treatment variable."),
                                call. = FALSE)
                    }
                    treat <- ptreat
                } else {
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
            } else if (length(propensity) == 2) {
                if (!hasArg(treat)) {
                    stop(gsub("\\s+", " ",
                              "Treatment variable is undetermined. Either
                           provide two-sided formula in the 'propensity'
                           argument, where the left hand variable is the
                           treatment variable, or declare the treatment variable
                           using the argument 'treat'."),
                         call. = FALSE)
                } else {
                    treat <- deparse(substitute(treat))
                }
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
            } else {
                stop(gsub("\\s+", " ",
                          "Treatment variable is undetermined. Either provide
                           two-sided formula in the 'propensity' argument,
                           where the left hand variable is the treatment
                           variable, or declare the treatment variable using
                           the argument 'treat'."),
                     call. = FALSE)
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

    ## Collect all variables declared, and remove unobserved variable
    ## from list
    allvars <- c(vars_y,
                 vars_formulas_x,
                 vars_formulas_z,
                 vars_subsets,
                 vars_mtr,
                 vars_weights,
                 vars_propensity)

    ## For the components, since they may be terms, we first collect
    ## all terms, and then break it down into variables.
    vars_components <- NULL
    if (userComponents) {
        for (comp in components) {

            compString <- try(gsub("\\s+", " ",
                                   Reduce(paste, deparse(comp))),
                              silent = TRUE)

            if (class(compString) != "try-error") {
                if (substr(compString, 1, 2) == "c(") {
                    vars_components <- c(vars_components,
                                         restring(comp,
                                                  substitute = FALSE))
                } else {
                    vars_components <- c(vars_components,
                                         restring(comp,
                                                  substitute = FALSE,
                                                  command = ""))
                }
            }
        }
    }
    vars_components <- vars_components[vars_components != '""']

    ## Break component terms into variables.
    vars_components_tmp <-
        paste(".qqq ~", paste(vars_components[vars_components != "components"],
                         collapse = " + "))
    if (! "intercept" %in% vars_components) {
        vars_components_tmp <- paste(vars_components_tmp, " - 1")
    }
    vars_components <- getXZ(as.formula(vars_components_tmp), components = TRUE)

    ## Collect all variables, and remove the variable name
    ## corresponding to the unobservable.
    allvars <- c(allvars, vars_components)
    allvars <- unique(allvars)
    allvars <- allvars[allvars != deparse(substitute(uname))]

    ## Fill in components list if necessary
    comp_filler <- lapply(terms_formulas_x,
                          function(x) as.character(unstring(x)))
    if (userComponents) {
        compMissing1 <- unlist(lapply(components, function(x) {
            Reduce(paste, deparse(x)) == ""
        }))
        compMissing2 <- unlist(lapply(components, function(x) x == ""))
        compMissing3 <- unlist(lapply(components, function(x) x == "c()"))
        compMissing <- as.logical(compMissing1 + compMissing2 + compMissing3)

        if (sum(compMissing) > 0) {
            components[compMissing] <- comp_filler[compMissing]
        }
    } else {
        components <- comp_filler
    }

    ## Check that all LATE variables are included in the propensity
    ## formula, if a propensity score formula is provided
    if (hasArg(target) && target == "late"  &&
        length(Formula::as.Formula(propensity))[1] != 0) {
        if (!all(late.Z %in% vars_propensity)) {
            stop(gsub("\\s+", " ",
                      "Variables in 'late.Z' argument must be contained
                       in the propensity score model."))
        }
        nLateX <- 0
        if (hasArg(late.X)) {
            nLateX <- length(late.X)
            if (!all(late.X %in% vars_propensity)) {
                stop(gsub("\\s+", " ",
                          "Variables in 'late.X' argument must be contained
                           in the propensity score model."))
            }
        }
        if (length(late.Z) + nLateX != length(vars_propensity) - 1) {
            stop(gsub("\\s+", " ",
                      "When estimating the LATE, all covariates in the
                       propensity model must be fixed using 'late.Z' and
                       'late.X'. Currently, not all variables are fixed."))
        }
    }

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

    ## Experimenting ------------------------------
    ## See what interacts with what
    print('original m0')
    print(origm0)
    print(splinesobj[[1]])
    
    stop('end of experimenting')
    ## End experiment -----------------------------
    
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
                                         "target.knots0", "target.knots1",
                                         "late.Z", "late.to", "late.from",
                                         "late.X", "eval.X"),
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

    if (hasArg(target) && target == "late") {
        estimateCall <- modcall(estimateCall,
                                newargs = list(late.Z = late.Z,
                                               late.to = late.to,
                                               late.from = late.from))
        if (hasArg(late.X)) {
            estimateCall <- modcall(estimateCall,
                                    newargs = list(late.X = late.X,
                                                   eval.X = eval.X))
        }
    }

    ## Estimate bounds
    if (point == FALSE) {

        origEstimate <- eval(estimateCall)

        ## Estimate bounds without resampling
        if (bootstraps == 0) {
            return(origEstimate)
        }
        ## Estimate bounds with resampling
        if (bootstraps > 0) {
            set.seed(seed)
            boundEstimates <- NULL

            b <- 1
            bootFailN <- 0
            bootFailNote <- ""
            bootFailIndex <- NULL

            if (!hasArg(bootstraps.m)) bootstraps.m <- nrow(data)
            if (bootstraps.m > nrow(data) && bootstraps.replace == FALSE) {
                stop(gsub("\\s+", " ",
                          "Cannot draw more observations than the number of rows
                           in the data set when 'bootstraps.replace = FALSE'."))
            }

            while (b <= bootstraps) {
                bootIDs  <- sample(seq(1, nrow(data)),
                                   size = bootstraps.m,
                                   replace = bootstraps.replace)
                bdata <- data[bootIDs, ]
                if (noisy == TRUE) {
                    message(paste0("Bootstrap iteration ", b, "..."))
                }
                bootCall <-
                    modcall(estimateCall,
                            dropargs = c("data", "noisy"),
                            newargs = list(data = quote(bdata),
                                           noisy = FALSE))
                bootEstimate <- try(eval(bootCall), silent = TRUE)
                if (is.list(bootEstimate)) {
                    boundEstimates  <- rbind(boundEstimates,
                                             bootEstimate$bound)
                    b <- b + 1
                    bootFailN <- 0
                    if (noisy == TRUE) {
                        message(paste0("    Audit count: ",
                                       bootEstimate$auditcount))
                        message(paste0("    Minimum criterion: ",
                                       fmtResult(bootEstimate$minobseq)))
                        message(paste0("    Bounds:",
                                       paste0("[",
                                              fmtResult(bootEstimate$bounds[1]),
                                              ", ",
                                              fmtResult(bootEstimate$bounds[2]),
                                              "]")))
                    }
                } else {
                    if (noisy == TRUE) {
                        stop(paste0("    Error:", bootEstimate))
                    }
                    bootFailN <- bootFailN + 1
                    bootFailIndex <- unique(c(bootFailIndex, b))
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

            if (ci.type == "both") {
                for (i in c("backward", "forward")) {
                    message(paste0("\nBootstrapped confidence intervals (",
                                   i, "):"))
                    for (j in 1:length(levels)) {
                        cistr <- paste0("[",
                                        fmtResult(ci[[i]][j, 1]),
                                        ", ",
                                        fmtResult(ci[[i]][j, 2]),
                                        "]")
                        message(paste0("    ",
                                       levels[j] * 100,
                                       "%: ",
                                       cistr))
                    }
                }
                message("\nBootstrapped p-values: ")
                message(paste0("    Backward: ", fmtResult(pvalue[1])))
                message(paste0("    Forward:  ", fmtResult(pvalue[2]), "\n"))
            } else {
                message(paste0("\nBootstrapped confidence intervals (",
                               ci.type, "):"))
                for (j in 1:length(levels)) {
                    cistr <- paste0("[",
                                    fmtResult(ci[j, 1]),
                                    ", ",
                                    fmtResult(ci[j, 2]),
                                    "]")
                    message(paste0("    ",
                                   levels[j] * 100,
                                   "%: ",
                                   cistr))
                }

                if (ci.type == "backward") {
                    message(paste0("\nBootstrapped p-value (backward): ",
                                   pvalue, "\n"))
                }
                if (ci.type == "forward") {
                    message(paste0("\nBootstrapped p-value (forward): ",
                                   pvalue, "\n"))
                }
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

    ## Point estimate without resampling
    if (point == TRUE & bootstraps == 0) {
        return(eval(estimateCall))
    }

    ## Point estimate with resampling
    if (point == TRUE & bootstraps > 0) {
        set.seed(seed)
        origEstimate <- eval(estimateCall)

        teEstimates  <- NULL
        mtrEstimates <- NULL
        propEstimates <- NULL

        b <- 1
        bootFailN <- 0
        bootFailNote <- ""
        bootFailIndex <- NULL

        if (!hasArg(bootstraps.m)) bootstraps.m <- nrow(data)
        if (bootstraps.m > nrow(data) && bootstraps.replace == FALSE) {
            stop(gsub("\\s+", " ",
                      "Cannot draw more observations than the number of rows
                           in the data set when 'bootstraps.replace = FALSE'."))
        }
        while (b <= bootstraps) {
            bootIDs  <- sample(seq(1, nrow(data)),
                                 size = bootstraps.m,
                                 replace = bootstraps.replace)
            bdata <- data[bootIDs, ]
            if (noisy == TRUE) {
                message(paste0("Bootstrap iteration ", b, "..."))
            }
            bootCall <-
                modcall(estimateCall,
                        dropargs = c("data", "noisy"),
                        newargs = list(data = quote(bdata),
                                       noisy = FALSE))
            bootEstimate <- try(eval(bootCall), silent = TRUE)
            if (is.list(bootEstimate)) {
                teEstimates  <- c(teEstimates, bootEstimate$pointestimate)
                mtrEstimates <- cbind(mtrEstimates, bootEstimate$mtr.coef)
                propEstimates <- cbind(propEstimates,
                                       bootEstimate$propensity$model$coef)
                b <- b + 1
                bootFailN <- 0
                if (noisy == TRUE) {
                    message(paste0("    Point estimate:",
                                   fmtResult(bootEstimate$pointestimate)))
                }
            } else {
                if (noisy == TRUE) {
                    message(paste0("    Error, resampling..."))
                }
                bootFailN <- bootFailN + 1
                bootFailIndex <- unique(c(bootFailIndex, b))
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

        ## Construct p-value
        pvalue <- (sum(teEstimates - origEstimate$pointestimate >=
                       abs(origEstimate$pointestimate)) +
                   sum(teEstimates - origEstimate$pointestimate <=
                       -abs(origEstimate$pointestimate))) / bootstraps

        ## Construct confidence intervals for various levels
        for (level in levels) {
            pLower <- (1 - level) / 2
            pUpper <- 1 - (1 - level) / 2
            probVec <- c(pLower, pUpper)

            ## Conf. int. 1: quantile method (same as percentile method)
            assign(paste0("ci1", level * 100),
                quantile(x = teEstimates,
                         probs = probVec,
                         type = 1))
            assign(paste0("mtrci1", level * 100),
                apply(mtrEstimates, 1, quantile,
                      probs = probVec,
                      type = 1))
            assign(paste0("propci1", level * 100),
                apply(propEstimates, 1, quantile,
                      probs = probVec,
                      type = 1))

            ## Conf. int. 2: percentile method using Z statistics
            tmpCi2 <- origEstimate$pointestimate +
                c(qnorm(pLower), qnorm(pUpper)) * bootSE
            names(tmpCi2) <- paste0(probVec * 100, "%")

            tmpMtrCi2 <- sweep(x = tcrossprod(c(qnorm(pLower),
                                                qnorm(pUpper)),
                                              mtrSE),
                               MARGIN = 2, origEstimate$mtr.coef, FUN = "+")
            tmpPropCi2 <- sweep(x = tcrossprod(c(qnorm(pLower),
                                                 qnorm(pUpper)),
                                               propSE), MARGIN = 2,
                                origEstimate$prop$model$coef, FUN = "+")
            colnames(tmpMtrCi2) <- colnames(get(paste0("mtrci1",
                                                       level * 100)))
            rownames(tmpMtrCi2) <- rownames(get(paste0("mtrci1",
                                                       level * 100)))
            colnames(tmpPropCi2) <- colnames(get(paste0("propci1",
                                                        level * 100)))
            rownames(tmpPropCi2) <- rownames(get(paste0("propci1",
                                                        level * 100)))

            assign(paste0("ci2", level * 100), tmpCi2)
            assign(paste0("mtrci2", level * 100), tmpMtrCi2)
            assign(paste0("propci2", level * 100), tmpPropCi2)
        }

        ## Prepare output
        output1 <- c(origEstimate,
                     list(pointestimate.se  = bootSE,
                          mtr.se = mtrSE,
                          prop.se = propSE,
                          pointestimate.bootstraps = teEstimates,
                          mtr.bootstraps = t(mtrEstimates)))
        output2 <- list()
        for (level in levels) {
            output2[[paste0("pointestimate.ci1.", level * 100)]] <-
                get(paste0("ci1", level * 100))
            output2[[paste0("mtr.ci1.", level * 100)]] <-
                t(get(paste0("mtrci1", level * 100)))
            output2[[paste0("prop.ci1.", level * 100)]] <-
                t(get(paste0("propci1", level * 100)))
            output2[[paste0("pointestimate.ci2.", level * 100)]] <-
                get(paste0("ci2", level * 100))
            output2[[paste0("mtr.ci2.", level * 100)]] <-
                t(get(paste0("mtrci2", level * 100)))
            output2[[paste0("prop.ci2.", level * 100)]] <-
                t(get(paste0("propci2", level * 100)))
        }
        output3 <- list(pvalue = pvalue,
                        bootstraps = bootstraps,
                        failed.bootstraps = length(bootFailIndex))
        output <- c(output1, output2, output3)

        message("\nBootstrapped confidence intervals (nonparametric):")
        for (level in levels) {
            ci1str <- get(paste0("ci1", level * 100))
            ci1str <- paste0("[",
                             fmtResult(ci1str[1]),
                             ", ",
                             fmtResult(ci1str[2]),
                             "]")
            message(paste0("    ",
                           level * 100,
                           "%: ",
                           ci1str))
        }
        message("\nBootstrapped confidence intervals (normal quantiles):")
        for (level in levels) {
            ci2str <- get(paste0("ci2", level * 100))
            ci2str <- paste0("[",
                             fmtResult(ci2str[1]),
                             ", ",
                             fmtResult(ci2str[2]),
                             "]")
            message(paste0("    ",
                           level * 100,
                           "%: ",
                           ci2str))
        }
        message(paste0("\nBootstrapped p-value: ",
                       fmtResult(pvalue), "\n"))
        return(output)
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
#'     Default is set to "logit".
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
#' @param vars_y character, variable name of observed outcome
#'     variable.
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
#' @param initgrid.nu number of evenly spread points in the interval
#'     [0, 1] of the unobservable u used to form the grid for imposing
#'     shape restrictions on the MTRs.
#' @param initgrid.nx number of evenly spread points of the covariates
#'     to use to form the grid for imposing shape restrictions on the
#'     MTRs.
#' @param audit.nx number of points on the covariates space to audit
#'     in each iteration of the audit procedure.
#' @param audit.nu number of points in the interval [0, 1],
#'     corresponding to the normalized value of the unobservable term,
#'     to audit in each iteration of the audit procedure
#' @param audit.add maximum number of points to add to the grids for
#'     imposing each kind of shape constraint. So if there are 5
#'     different kinds of shape constraints, there can be at most
#'     \code{audit.add * 5} additional points added to the grid.
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
#' @param point.eyeweight boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} if GMM point estimate should use the identity
#'     weighting matrix (i.e. one-step GMM).
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @param seed integer, the seed that determines the random grid in
#'     the audit procedure.
#' @return Returns a list of results from throughout the estimation
#'     procedure. This includes all IV-like estimands; the propensity
#'     score model; bounds on the treatment effect; the estimated
#'     expectations of each term in the MTRs; the components and
#'     results of the LP problem.
ivmteEstimate <- function(ivlike, data, subset, components,
                          propensity, link = "logit", treat, m0, m1,
                          vars_y, vars_mtr, terms_mtr0, terms_mtr1,
                          splinesobj, uname = u, target,
                          target.weight0, target.weight1,
                          target.knots0 = NULL, target.knots1 = NULL,
                          late.Z, late.from, late.to, late.X, eval.X,
                          genlate.lb, genlate.ub, obseq.tol = 0.05,
                          initgrid.nu = 10, initgrid.nx = 20, audit.nx = 2500,
                          audit.nu = 25, audit.add = 100, audit.max = 25,
                          audit.tol = 1e-08, m1.ub, m0.ub, m1.lb,
                          m0.lb, mte.ub, mte.lb, m0.dec, m0.inc,
                          m1.dec, m1.inc, mte.dec, mte.inc,
                          lpsolver = NULL, point = FALSE,
                          point.eyeweight = FALSE,
                          noisy = TRUE, seed = 12345) {

    call <- match.call(expand.dots = FALSE)

    if (classFormula(ivlike)) ivlike <- c(ivlike)

    ## Character arguments will be converted to lowercase
    if (hasArg(lpsolver)) lpsolver <- tolower(lpsolver)
    if (hasArg(target))   target   <- tolower(target)
    if (hasArg(link))     link     <- tolower(link)
    if (hasArg(ci.type))  ci.type  <- tolower(ci.type)


    ##---------------------------
    ## 1. Obtain propensity scores
    ##---------------------------

    if (noisy == TRUE) {
        message("Obtaining propensity scores...")
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
        message("Generating target moments...")
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

    if (classList(ivlike)) {
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
                                  noisy = noisy,
                                  ivn = ivlikeCounter)
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
                                  noisy = noisy,
                                  ivn = ivlikeCounter)
            }
            ivlikeCounter <- ivlikeCounter + 1

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
                                 identity = point.eyeweight,
                                 noisy = noisy)

        return(list(sset  = sset,
                    gstar = list(g0 = colMeans(gstar0),
                                 g1 = colMeans(gstar1)),
                    propensity = pmodel,
                    pointestimate = gmmResult$pointestimate,
                    Jtest = gmmResult$Jtest,
                    bounds = c(gmmResult$pointestimate,
                               gmmResult$pointestimate),
                    mtr.coef = gmmResult$coef))
    }

    ##---------------------------
    ## 4. Define constraint matrices using the audit
    ##---------------------------

    if (noisy == TRUE) {
        message("Performing audit procedure...")
    }

    audit.args <- c("uname", "initgrid.nu", "initgrid.nx",
                    "audit.nx", "audit.nu", "audit.add",
                    "audit.max", "audit.tol",
                    "m1.ub", "m0.ub",
                    "m1.lb", "m0.lb",
                    "mte.ub", "mte.lb", "m0.dec",
                    "m0.inc", "m1.dec", "m1.inc", "mte.dec",
                    "mte.inc", "obseq.tol", "noisy", "seed")

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
    if (noisy) {
        message(paste0("Bounds on the target parameter: [",
                       fmtResult(audit$min), ", ", fmtResult(audit$max), "]\n"))
    }
    ## include additional output material
    return(list(sset  = sset,
                gstar = list(g0 = gstar0,
                             g1 = gstar1,
                             w0 = targetGammas$w0,
                             w1 = targetGammas$w1),
                propensity = pmodel,
                pointestimate = NULL,
                bounds = c(audit$min, audit$max),
                lpresult =  audit$lpresult,
                auditgrid = audit$gridobj,
                auditcount = audit$auditcount,
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

    if (hasArg(target)) target   <- tolower(target)

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
            w1 <- wlate1(data, late.from, late.to, late.Z,
                         pmodobj$model, late.X, eval.X)
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
                message("    Integrating terms for control group...")
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
                message("    Integrating terms for treated group...")
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
                            "    Integrating non-spline terms for control group...")
                    } else {
                        message(
                            "    Integrating non-spline terms for treated group...")
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
    output <- list(gstar0 = gstar0,
                   gstar1 = gstar1,
                   xindex0 = xindex0,
                   xindex1 = xindex1,
                   uexporder0 = uexporder0,
                   uexporder1 = uexporder1)
    if (hasArg(target)) {
        output$w1 <- w1
        output$w0 <- w0
    }
    return(output)
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
#' @param ivn integer, the number indicating which IV specification
#'     the component corresponds to.
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
                        yvar, dvar, noisy = TRUE, ivn = NULL) {
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
        sweight0 <- list(lb = pmodobj,
                         ub = 1,
                         multiplier = sest$sw0[, j])
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
        sweight1 <- list(lb = 0,
                         ub = pmodobj,
                         multiplier = sest$sw1[, j])
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
            sset[[paste0("s", scount)]] <- list(ivspec = ivn,
                                                beta = sest$beta[j],
                                                g0 = c(gs0, gsSpline0),
                                                g1 = c(gs1, gsSpline1),
                                                w0 = sweight0,
                                                w1 = sweight1)
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
            sset[[paste0("s", scount)]] <- list(ivspec = ivn,
                                                beta = sest$beta[j],
                                                g0 = cbind(gs0, gsSpline0),
                                                g1 = cbind(gs1, gsSpline1),
                                                ys = yvec,
                                                w0 = sweight0,
                                                w1 = sweight1)
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
#' condition is then set up as a two-step GMM problem. Solving this
#' GMM problem recovers the coefficients on the MTR functions m0 and
#' m1. Combining these coefficients with the target gamma moments
#' allows us to estimate the target treatment effect.
#' @param sset a list of lists constructed from the function
#'     \link{genSSet}. Each inner list should include a coefficient
#'     corresponding to a term in an IV specification, a matrix of the
#'     estimates of the gamma moments conditional on (X, Z) for d = 0,
#'     and a matrix of the estimates of the gamma moments conditional
#'     on (X, Z) for d = 1. The column means of the last two matrices
#'     is what is used to generate the gamma moments.
#' @param gstar0 vector, the target gamma moments for d = 0.
#' @param gstar1 vector, the target gamma moments for d = 1.
#' @param identity boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} if GMM point estimate should use the identity
#'     weighting matrix (i.e. one-step GMM).
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
#' ## Obtain point estimates using GMM
#' gmmEstimate(sset = sSet$sset,
#'             gstar0 = targetGamma$gstar0,
#'             gstar1 = targetGamma$gstar1)
#'
#' @export
gmmEstimate <- function(sset, gstar0, gstar1, identity = FALSE, noisy = TRUE) {
    gn0 <- ncol(gstar0)
    gn1 <- ncol(gstar1)
    if ((gn0 + gn1) > length(sset)) {
        stop(gsub("\\s+", " ",
                  paste0("System is underidentified: there are ",
                         (gn0 + gn1),
                         " unknown MTR coefficients and ",
                         length(sset),
                         " moment conditions defined by IV-like
                         specifications. Either expand the number of
                         IV-like specifications, or modify m0 and m1.")))
    }

    ## This function constructs the matrix to be fed into the GMM
    ## estimator to construct the moment conditions.
    momentMatrix <- function() {
        momentMatrix <- NULL
        momentNames <- NULL
        for (s in 1:length(sset)) {
            momentMatrix <- cbind(momentMatrix,
                                  sset[[s]]$ys , sset[[s]]$g0,  sset[[s]]$g1)
            momentNames <- c(momentNames,
                             paste0("s", s, "y"),
                             paste0("s", s, "g0", seq(1, gn0)),
                             paste0("s", s, "g1", seq(1, gn1)))
        }
        colnames(momentMatrix) <- momentNames
        return(momentMatrix)
    }

    ## This function defines the moment conditions for the GMM
    ## estimator. The argument 'theta' is for a vector storing the
    ## parameters of interest in the following order: (i) coefficients
    ## for m0; (ii) coefficients for m1.
    momentConditions <- function(theta, data) {
        momentMatrix <- NULL
        for (s in 1:length(sset)) {
            yvec <- data[, paste0("s", s, "y")]
            gmat <- data[, c( paste0("s", s, "g0", seq(1, gn0)),
                             paste0("s", s, "g1", seq(1, gn1)))]
            momentMatrix <- cbind(momentMatrix,
                                  yvec - gmat %*% theta[1:(gn0 + gn1)])
        }
        return(momentMatrix)
    }

    ## Perform GMM
    if (identity == FALSE) {
        gmmObj <- gmm::gmm(momentConditions, x = momentMatrix(),
                           t0 = rep(0, times = (gn0 + gn1)),
                           prewhite = 1)
    } else {
        gmmObj <- gmm::gmm(momentConditions, x = momentMatrix(),
                           t0 = rep(0, times = (gn0 + gn1)),
                           prewhite = 1, wmatrix = "ident")
    }
    theta <- gmmObj$coefficients
    if (length(sset) > gn0 + gn1) {
        Jtest <- gmmObj$objective * nrow(gstar0)
        Jtest <- c(Jtest,
                   1 - pchisq(Jtest, df = length(sset) - gn0 - gn1))
        names(Jtest) <- c("J-statistic", "p-value")
    } else {
        Jtest <- NULL
    }

    ## Construct point estimate and CI of TE
    names(theta) <- c(paste0("m0.", colnames(gstar0)),
                      paste0("m1.", colnames(gstar1)))
    pointestimate <- sum(c(colMeans(gstar0), colMeans(gstar1)) * theta)
    if (noisy == TRUE) {
        message()
        message(paste0("Point estimate of the target parameter: ",
                       round(pointestimate, 4), "\n"))
    }
    return(list(pointestimate = as.numeric(pointestimate),
                coef = theta,
                Jtest = Jtest))
}

#' Format result for display
#'
#' This function simply takes a number and formats it for being
#' displayed. Numbers less than 1 in absolute value are rounded to 6
#' significant figure. Numbers larger than
#'
#' @param x The scalar to be formated
#' @return A scalar.
fmtResult <- function(x) {
    if (abs(x) < 1) {
        fx <- signif(x, digits = 7)
    } else if (abs(x) >= 1 & abs(x) < 1e+7) {
        fx <- signif(round(x, digits = 4), digits = 7)
    } else {
        fx <- formatC(x, format = "e", digits = 7)
    }
    if (is.numeric(fx)) fx <- as.character(fx)
    return(fx)
}
