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
#'     function for control group.
#' @param m1 one-sided formula for marginal treatment response
#'     function for treated group.
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
#' @param grid.Nu number of evenly spread points in the interval [0,
#'     1] of the unobservable u used to form the grid for imposing
#'     shape restrictions on the MTRs.
#' @param grid.Nx number of evenly spread points of the covariates to
#'     use to form the grid for imposing shape restrictions on the
#'     MTRs.
#' @param audit.Nx number of points on the covariates space to audit
#'     in each iteration of the audit procedure.
#' @param audit.Nu number of points in the interval [0, 1],
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
ivmte <- function(ivlike, data, subset, components, propensity, link,
                  treat, m0, m1, uname = u, target, target.weight0,
                  target.weight1, target.knots0 = NULL, target.knots1 = NULL,
                  late.Z, late.from, late.to, late.X,
                  eval.X, genlate.lb, genlate.ub, obseq.tol = 1,
                  grid.Nu = 20, grid.Nx = 20, audit.Nx = 2,
                  audit.Nu = 3, audit.max = 5, audit.tol = 1e-08, m1.ub,
                  m0.ub, m1.lb, m0.lb, mte.ub, mte.lb, m0.dec, m0.inc,
                  m1.dec, m1.inc, mte.dec, mte.inc, lpsolver = NULL) {
    
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

        if (hasArg(components) & !is.null(components)) {
            userComponents <- TRUE
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
            userComponents <- FALSE
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

            ## Conver all entries into lists
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
    
    ## Audit checks
    if (!(is.numeric(obseq.tol) & obseq.tol >= 0)) {
        stop("Cannot set obseq.tol below 0.")
    }
    if (!((grid.Nu %% 1 == 0) & grid.Nu >= 2)) {
        stop("grid.Nu must be an integer greater than or equal to 2.")
    }
    if (!((grid.Nx %% 1 == 0) & grid.Nx >= 0)) {
        stop("grid.Nx must be an integer greater than or equal to 0.")
    }

    if (!((audit.Nx %% 1 == 0) & audit.Nx > 0)) {
        stop("audit.Nx must be an integer greater than or equal to 1.")
    }

    if (hasArg(m0.dec) | hasArg(m0.inc) |
        hasArg(m1.dec) | hasArg(m1.inc) |
        hasArg(mte.dec) | hasArg(mte.inc)) {
        if (!((audit.Nu %% 1 == 0) & audit.Nu > 0) | audit.Nu < 2) {
            stop("audit.Nu must be an integer greater than or equal to 2.")
        }
    } else {
        if (!((audit.Nu %% 1 == 0) & audit.Nu > 0)) {
            stop("audit.Nu must be an integer greater than or equal to 1.")
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
            treat <- all.vars(propensity)[1]
            vars_propensity <- all.vars(propensity)

            if (hasArg(treat)) {
                if (treat != deparse(substitute(treat))) {
                    warning(gsub("\\s+", " ",
                                 "Variable listed in 'treat' argument differs
                                 from dependent variable in propensity score
                                 formula. Dependent variable from propensity
                                 score formula will be used as the treatment
                                 variable."),
                            call. = FALSE)
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
    cdata <- data

    ##---------------------------
    ## 2. Obtain propensity scores
    ##---------------------------

    message("Obtaining propensity scores...\n")

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

    message("Generating target moments...\n")

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
            w1 <- wate1.mst(data)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "att") {
            w1 <- watt1.mst(data, mean(data[[treat]]), pmodel$phat)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "atu") {
            w1 <- watu1.mst(data, 1- mean(data[[treat]]), pmodel$phat)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "late") {
            if (!hasArg(late.X)) {
                late.X <- NULL
                eval.X <- NULL
            }
            w1 <- wlate1.mst(data, late.from, late.to, substitute(late.Z),
                             pmodel$model, substitute(late.X), eval.X)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "genlate") {
            w1 <- wgenlate1.mst(data, genlate.lb, genlate.ub)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else {
            stop("Unrecognized target parameter.")
        }

        ## Integrate m0 and m1 functions
        if (!is.null(m0)) {
            message("    Integrating terms for control group...\n")
            pm0 <- eval(as.call(m0call))
            gstar0 <- gengamma.mst(pm0, w0$lb, w0$ub, w0$mp)
        } else {
            gstar0 <- NULL
            pm0 <- NULL
        }

        if (!is.null(m1)) {
            message("    Integrating terms for treated group...\n")
            pm1 <- eval(as.call(m1call))
            gstar1 <- gengamma.mst(pm1, w1$lb, w1$ub, w1$mp)
        } else {
            gstar1 <- NULL
            pm1 <- NULL
        }

        gstarSpline0 <- genGammaSplines.mst(splines = splinesobj[[1]],
                                            data = cdata,
                                            lb = w0$lb,
                                            ub = w0$ub,
                                            multiplier = w0$mp,
                                            d = 0)

        gstarSpline1 <- genGammaSplines.mst(splines = splinesobj[[2]],
                                            data = cdata,
                                            lb = w1$lb,
                                            ub = w1$ub,
                                            multiplier = w1$mp,
                                            d = 1)
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
                if (d == 0) {
                    message(
                      "    Integrating non-spline terms for control group...\n")
                } else {
                    message(
                      "    Integrating non-spline terms for treated group...\n")
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

                    gamma <- gamma + gengamma.mst(pm,
                                                  lb = lb,
                                                  ub = ub,
                                                  multiplier = weights)
                }

                assign(paste0("gstar", d), gamma)
                assign(paste0("pm", d), pm)
            } else {
                assign(paste0("gstar", d), NULL)
                assign(paste0("pm", d), NULL)
            }

            ## Integrate splines terms

            if (d == 0) {
                message(
                    "    Integrating spline terms for control group...\n")
            } else {
                message(
                    "    Integrating spline terms for treated group...\n")
            }

            noSplineMtr <- splinesobj[[d + 1]]

            if (!is.null(noSplineMtr$splineslist)) {

                basisLengths <- sapply(names(noSplineMtr$splineslist),
                                       function(x) {
                    argList <- eval(parse(text = gsub("uSplines\\(",
                                                      "list(", x)))
                    basisLength <- length(argList$knots) +  argList$degree
                    if(isTRUE(argList$intercept == TRUE) |
                       is.null(argList$intercept)) basisLength <-
                                                       basisLength + 1
                    return(basisLength)
                })

                gammaSplines <- rep(0, sum(basisLengths))

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

                    gammaSplines <- gammaSplines +
                        genGammaSplines.mst(splines = noSplineMtr,
                                            data = cdata,
                                            lb = lb,
                                            ub = ub,
                                            multiplier = weights,
                                            d = d)
                }
                assign(paste0("gstarSpline", d), gammaSplines)
            } else {
                assign(paste0("gstarSpline", d), NULL)
            }
        }
    }

    gstar0 <- c(gstar0, gstarSpline0)
    gstar1 <- c(gstar1, gstarSpline1)

    ##---------------------------
    ## 4. Generate moments/gamma terms for IV-like estimands
    ##---------------------------

    message("Generating IV-like moments...")

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
        setobj <- gensset.mst(data = cdata,
                              sset = sset,
                              sest = sest,
                              splinesobj = splinesobj,
                              pmodobj = pmodel$phat[subset_index],
                              pm0 = pm0,
                              pm1 = pm1,
                              ncomponents = ncomponents,
                              scount = scount,
                              subset_index = subset_index)

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
            setobj <- gensset.mst(data = cdata,
                                  sset = sset,
                                  sest = sest,
                                  splinesobj = splinesobj,
                                  pmodobj = pmodobj,
                                  pm0 = pm0,
                                  pm1 = pm1,
                                  ncomponents = ncomponents,
                                  scount = scount,
                                  subset_index = subset_index)

            ## Update set of moments (gammas)
            sset <- setobj$sset
            scount <- setobj$scount
        }
    } else {
        stop(gsub("\\s+", " ",
                  "'ivlike' argument must either be a formula or a vector of
                  formulas."))
    }

    ##---------------------------
    ## 5. Define constraint matrices using the audit
    ##---------------------------

    if (obseq.tol > 0) {
        message("\nPerforming audit procedure...\n")
    }

    audit.args <- c("uname", "grid.Nu", "grid.Nx",
                    "audit.Nx", "audit.Nu", "audit.max", "audit.tol",
                    "m1.ub", "m0.ub",
                    "m1.lb", "m0.lb", "mte.ub", "mte.lb", "m0.dec",
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
                                  newargs = list(m1.ub = maxy))
        }
        if (!hasArg(m0.ub)) {
            audit_call <- modcall(audit_call,
                                  newargs = list(m0.ub = maxy))
        }
    }
    if (!hasArg(m1.lb) | !hasArg(m0.lb)) {
        miny <- min(cdata[, vars_y])
        if (!hasArg(m1.lb)) {
            audit_call <- modcall(audit_call,
                                  newargs = list(m1.lb = miny))
        }
        if (!hasArg(m0.lb)) {
            audit_call <- modcall(audit_call,
                                  newargs = list(m0.lb = miny))
        }
    }

    ##---------------------------
    ## 6. Obtain the bounds
    ##---------------------------

    audit <- eval(audit_call)
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
                minobseq = audit$minobseq))
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
#' @return A list containing the point estimate for the IV regression,
#'     and the expectation of each monomial term in the MTR.
gensset.mst <- function(data, sset, sest, splinesobj, pmodobj, pm0, pm1,
                        ncomponents, scount, subset_index) {

    for (j in 1:ncomponents) {
        message(paste0("    Moment ", scount, "..."))

        if (!is.null(pm0)) {
            gs0 <- gengamma.mst(monomials = pm0,
                                lb = pmodobj,
                                ub = 1,
                                multiplier = sest$sw0[, j],
                                subset = subset_index)
        } else {
            gs0 <- NULL
        }

        if (!is.null(pm1)) {
            gs1 <- gengamma.mst(monomials = pm1,
                                lb = 0,
                                ub = pmodobj,
                                multiplier = sest$sw1[, j],
                                subset = subset_index)
        } else {
            gs1 <- NULL
        }

        gsSpline0 <- genGammaSplines.mst(splines = splinesobj[[1]],
                                         data = data,
                                         lb = pmodobj,
                                         ub = 1,
                                         multiplier = sest$sw0[, j],
                                         subset = subset_index,
                                         d = 0)

        gsSpline1 <- genGammaSplines.mst(splines = splinesobj[[2]],
                                         data = data,
                                         lb = 0,
                                         ub = pmodobj,
                                         multiplier = sest$sw1[, j],
                                         subset = subset_index,
                                         d = 1)

        ## generate components of constraints
        sset[[paste0("s", scount)]] <- list(beta = sest$beta[j],
                                            g0 = c(gs0, gsSpline0),
                                            g1 = c(gs1, gsSpline1))

        ## update counter (note scount is not referring
        ## to the list of IV regressions, but the components
        ## from the IV regressions)
        scount <- scount + 1
    }

    return(list(sset = sset, scount = scount))
}

