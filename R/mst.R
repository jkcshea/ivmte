utils::globalVariables("u")
#' Instrumental Variables: Extrapolation by Marginal Treatment Effects
#'
#' This function provides a general framework for using the marginal
#' treatment effect (MTE) to extrapolate. The model is the same binary
#' treatment instrumental variable (IV) model considered by
#' \href{https://doi.org/10.2307/2951620}{Imbens and Angrist (1994)}
#' and \href{https://doi.org/10.1111/j.1468-0262.2005.00594.x}{Heckman
#' and Vytlacil (2005)}. The framework on which this function is based
#' was developed by \href{https://doi.org/10.3982/ECTA15463}{Mogstad,
#' Santos and Torgovitsky (2018)}. See also the recent survey paper on
#' extrapolation in IV models by
#' \href{https://doi.org/10.1146/annurev-economics-101617-041813}{Mogstad
#' and Torgovitsky (2018)}. A detailed description of the module and
#' its features can be found in
#' \href{https://a-torgovitsky.github.io/shea-torgovitsky.pdf}{Shea
#' and Torgovitsky (2019)}.
#'
#' @import methods stats utils
#'
#' @param data \code{data.frame} or \code{data.table} used to estimate
#'     the treatment effects.
#' @param target character, target parameter to be
#'     estimated. Currently function allows for ATE (\code{'ate'}),
#'     ATT (\code{'att'}), ATU (\code{'atu'}), LATE (\code{'late'}),
#'     and generalized LATE (\code{'genlate'}).
#' @param late.from a named vector, or a list, declaring the baseline
#'     set of values of Z used to define the LATE. The name associated
#'     with each value should be the name of the corresponding
#'     variable.
#' @param late.to a named vector, or a list, declaring the comparison
#'     set of values of Z used to define the LATE. The name associated
#'     with each value should be the name of the corresponding
#'     variable.
#' @param late.X a named vector, or a list, declaring the values at
#'     which to condition on. The name associated with each value
#'     should be the name of the corresponding variable.
#' @param genlate.lb lower bound value of unobservable \code{u} for
#'     estimating the generalized LATE.
#' @param genlate.ub upper bound value of unobservable \code{u} for
#'     estimating the generalized LATE.
#' @param target.weight0 user-defined weight function for the control
#'     group defining the target parameter. A list of functions can be
#'     submitted if the weighting function is in fact a spline. The
#'     arguments of the function should be variable names in
#'     \code{data}. If the weight is constant across all observations,
#'     then the user can instead submit the value of the weight
#'     instead of a function.
#' @param target.weight1 user-defined weight function for the treated
#'     group defining the target parameter. See \code{target.weight0}
#'     for details.
#' @param target.knots0 user-defined set of functions defining the
#'     knots associated with spline weights for the control group. The
#'     arguments of the function should consist only of variable names
#'     in \code{data}. If the knots are constant across all
#'     observations, then the user can instead submit the vector of
#'     knots instead of a function.
#' @param target.knots1 user-defined set of functions defining the
#'     knots associated with spline weights for the treated group. See
#'     \code{target.knots0} for details.
#' @param m0 one-sided formula for the marginal treatment response
#'     function for the control group. Splines may also be
#'     incorporated using the expression \code{uSpline}, e.g.
#'     \code{uSpline(degree = 2, knots = c(0.4, 0.8), intercept =
#'     TRUE)}. The \code{intercept} argument may be omitted, and is
#'     set to \code{TRUE} by default.
#' @param m1 one-sided formula for marginal treatment response
#'     function for treated group. Splines can also be incorporated
#'     using the expression "uSplines(degree, knots, intercept)". The
#'     \code{intercept} argument may be omitted, and is set to
#'     \code{TRUE} by default.
#' @param uname variable name for the unobservable used in declaring
#'     the MTRs. The name can be provided with or without quotation
#'     marks.
#' @param m1.ub numeric value for upper bound on MTR for the treated
#'     group. By default, this will be set to the largest value of the
#'     observed outcome in the estimation sample.
#' @param m0.ub numeric value for upper bound on MTR for the control
#'     group. By default, this will be set to the largest value of the
#'     observed outcome in the estimation sample.
#' @param m1.lb numeric value for lower bound on MTR for the treated
#'     group. By default, this will be set to the smallest value of
#'     the observed outcome in the estimation sample.
#' @param m0.lb numeric value for lower bound on MTR for the control
#'     group. By default, this will be set to the smallest value of
#'     the observed outcome in the estimation sample.
#' @param mte.ub numeric value for upper bound on treatment effect
#'     parameter of interest.
#' @param mte.lb numeric value for lower bound on treatment effect
#'     parameter of interest.
#' @param m0.dec logical, set to \code{FALSE} by default. Set equal to
#'     \code{TRUE} if the MTR for the control group should be weakly
#'     monotone decreasing.
#' @param m0.inc logical, set to \code{FALSE} by default. Set equal to
#'     \code{TRUE} if the MTR for the control group should be weakly
#'     monotone increasing.
#' @param m1.dec logical, set to \code{FALSE} by default. Set equal to
#'     \code{TRUE} if the MTR for the treated group should be weakly
#'     monotone decreasing.
#' @param m1.inc logical, set to \code{FALSE} by default. Set equal to
#'     \code{TRUE} if the MTR for the treated group should be weakly
#'     monotone increasing.
#' @param mte.dec logical, set to \code{FALSE} by default. Set equal
#'     to \code{TRUE} if the MTE should be weakly monotone decreasing.
#' @param mte.inc logical, set to \code{FALSE} by default. Set equal
#'     to \code{TRUE} if the MTE should be weakly monotone increasing.
#' @param ivlike formula or vector of formulas specifying the
#'     regressions for the IV-like estimands. Which coefficients to
#'     use to define the constraints determining the treatment effect
#'     bounds (alternatively, the moments determining the treatment
#'     effect point estimate) can be selected in the argument
#'     \code{components}.
#' @param components a list of vectors of the terms in the regression
#'     specifications to include in the set of IV-like estimands. No
#'     terms should be in quotes. To select the intercept term,
#'     include the name \code{intercept}. If the factorized
#'     counterpart of a variable is included in the IV-like
#'     specifications, e.g. \code{factor(x)} where \code{x = 1, 2, 3},
#'     the user can select the coefficients for specific factors by
#'     declaring the components \code{factor(x)-1, factor(x)-2,
#'     factor(x)-3}. See \code{\link{l}} on how to input the
#'     argument. If no components for a IV specification are given,
#'     then all coefficients from that IV specification will be used
#'     to define constraints in the partially identified case, or to
#'     define moments in the point identified case.
#' @param subset a single subset condition or list of subset
#'     conditions corresponding to each regression specified in
#'     \code{ivlike}. The input must be logical. See \code{\link{l}}
#'     on how to input the argument. If the user wishes to select
#'     specific rows, construct a binary variable in the data set, and
#'     set the condition to use only those observations for which the
#'     binary variable is 1, e.g. the binary variable is \code{use},
#'     and the subset condition is \code{use == 1}.
#' @param propensity formula or variable name corresponding to
#'     propensity to take up treatment. If a formula is declared, then
#'     the function estimates the propensity score according to the
#'     formula and link specified in \code{link}. If a variable name
#'     is declared, then the corresponding column in the data is taken
#'     as the vector of propensity scores. A variable name can be
#'     passed either as a string (e.g \code{propensity = 'p'}). , a
#'     variable (e.g. \code{propensity = p}), or a one-sided formula
#'     (e.g. \code{propensity = ~p}.
#' @param link character, name of link function to estimate propensity
#'     score. Can be chosen from \code{'linear'}, \code{'probit'}, or
#'     \code{'logit'}. Default is set to \code{'logit'}.
#' @param treat variable name for treatment indicator. The name can be
#'     provided with or without quotation marks.
#' @param lpsolver character, name of the linear programming package
#'     in R used to obtain the bounds on the treatment effect. The
#'     function supports \code{'gurobi'}, \code{'cplexapi'},
#'     \code{'lpsolveapi'}.
#' @param lpsolver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @param lpsolver.presolve boolean, default set to \code{TRUE}. Set
#'     this parameter to \code{FALSE} if presolve should be turned off
#'     for the LP problems.
#' @param lpsolver.options.criterion list, each item of the list
#'     should correspond to an option specific to the LP solver
#'     selected. These options are specific for finding the minimum
#'     criterion.
#' @param lpsolver.options.bounds list, each item of the list should
#'     correspond to an option specific to the LP solver
#'     selected. These options are specific for finding the bounds.
#' @param criterion.tol tolerance for violation of observational
#'     equivalence, set to 0 by default. Statistical noise may
#'     prohibit the theoretical LP problem from being feasible. That
#'     is, there may not exist a set of coefficients on the MTR that
#'     are observationally equivalent with regard to the IV-like
#'     regression coefficients. The function therefore first estimates
#'     the minimum violation of observational equivalence. This is
#'     reported in the output under the name 'minimum criterion'. The
#'     constraints in the LP problem pertaining to observational
#'     equivalence are then relaxed by the amount \code{minimum
#'     criterion * (1 + criterion.tol)}. Set \code{criterion.tol} to a
#'     value greater than 0 to allow for more conservative bounds.
#' @param initgrid.nx integer determining the number of points of the
#'     covariates used to form the initial constraint grid for
#'     imposing shape restrictions on the MTRs.
#' @param initgrid.nu integer determining the number of points in the
#'     open interval (0, 1) drawn from a Halton sequence. The end
#'     points 0 and 1 are additionally included. These points are used
#'     to form the initial constraint grid for imposing shape
#'     restrictions on the \code{u} components of the MTRs.
#' @param audit.nx integer determining the number of points on the
#'     covariates space to audit in each iteration of the audit
#'     procedure.
#' @param audit.nu integer determining the number of points in the
#'     interval [0, 1], corresponding to space of unobservable
#'     \code{u}, to audit in each iteration of the audit procedure.
#' @param audit.add maximum number of points to add to the initial
#'     constraint grid for imposing each kind of shape constraint. For
#'     example, if there are 5 different kinds of shape constraints,
#'     there can be at most \code{audit.add * 5} additional points
#'     added to the constraint grid.
#' @param audit.max maximum number of iterations in the audit
#'     procedure.
#' @param audit.tol feasibility tolerance when performing the
#'     audit. By default to set to be equal to the Gurobi
#'     (\code{lpsolver = "gurobi"}) and CPLEX (\code{lpsolver =
#'     "cplexapi"}) feasibility toleraence, which is set to
#'     \code{1e-06} by default.  If the LP solver is lp_solve
#'     (\code{lpsolver = "lpsolveapi"}), this parameter is set to
#'     \code{1e-06} by default. This parameter should only be changed
#'     if the feasibility tolerance of the LP solver is changed, or if
#'     numerical issues result in discrepancies between the LP
#'     solver's feasibility check and the audit.
#' @param point boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} if it is believed that the treatment effects are
#'     point identified. If set to \code{TRUE}, then a two-step GMM
#'     procedure is implemented to estimate the treatment
#'     effects. Shape constraints on the MTRs will be ignored under
#'     point identification.
#' @param point.eyeweight boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} if the GMM point estimate should use the identity
#'     weighting matrix (i.e. one-step GMM).
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
#' @param levels vector of real numbers between 0 and 1. Values
#'     correspond to the level of the confidence intervals constructed
#'     via bootstrap.
#' @param ci.type character, default set to \code{'both'}. Set to
#'     \code{'forward'} to construct the forward confidence interval
#'     for the treatment effect bound. Set to \code{'backward'} to
#'     construct the backward confidence interval for the treatment
#'     effect bound. Set to \code{'both'} to construct both types of
#'     confidence intervals.
#' @param specification.test boolean, default set to
#'     \code{TRUE}. Function performs a specificaiton test for the
#'     partially identified case when \code{bootstraps > 0}.
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @param smallreturnlist boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} to exclude large intermediary components
#'     (i.e. propensity score model, LP model, bootstrap iterations)
#'     from being included in the return list.
#' @param seed integer, the seed that determines the random grid in
#'     the audit procedure.
#' @param debug boolean, indicates whether or not the function should
#'     provide output when obtaining bounds. The option is only
#'     applied when \code{lpsolver = 'gurobi'}. The output provided is
#'     the same as what the Gurobi API would send to the console.
#' @return Returns a list of results from throughout the estimation
#'     procedure. This includes all IV-like estimands; the propensity
#'     score model; bounds on the treatment effect; the estimated
#'     expectations of each term in the MTRs; the components and
#'     results of the LP problem.
#'
#' @details The return list includes the following objects.
#'     \describe{
#' \item{sset}{a list of all the coefficient estimates and weights
#'     corresponding to each element in the S-set.}
#' \item{gstar}{a list containing the estimate of the weighted means
#' for each component in the MTRs. The weights are determined by the
#' target parameter declared in \code{target}, or the weights defined
#' by \code{target.weight1}, \code{target.knots1},
#' \code{target.weight0}, \code{target.knots0}.}
#' \item{gstar.weights}{a list containing the target weights used to
#' estimate \code{gstar}.}
#' \item{gstar.coef}{a list containing the coefficients on the treated
#' and control group MTRs.}
#' \item{propensity}{the propensity score model. If a variable is fed
#' to the \code{propensity} argument when calling \code{ivmte}, then
#' the returned object is a list containing the name of variable given by the
#' user, and the values of that variable used in estimation.}
#' \item{bounds}{a vector with the estimated lower and upper bounds of
#' the target treatment effect.}
#' \item{lpresult}{a list containing the LP model, and the full output
#' from solving the LP problem.}
#' \item{audit.grid}{the audit grid on which all shape constraints
#' were satisfied.}
#' \item{audit.count}{the number of audits required until there were
#' no more violations.}
#' \item{audit.criterion}{the minimum criterion.}
#' \item{splinesdict}{a list including the specifications of each
#' spline declared in each MTR.}
#' }
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#'
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
ivmte <- function(data, target, late.from, late.to, late.X,
                  genlate.lb, genlate.ub, target.weight0 = NULL,
                  target.weight1 = NULL, target.knots0 = NULL,
                  target.knots1 = NULL, m0, m1, uname = u, m1.ub,
                  m0.ub, m1.lb, m0.lb, mte.ub, mte.lb, m0.dec, m0.inc,
                  m1.dec, m1.inc, mte.dec, mte.inc, ivlike,
                  components, subset, propensity, link = 'logit',
                  treat, lpsolver = NULL, lpsolver.options,
                  lpsolver.presolve,
                  lpsolver.options.criterion, lpsolver.options.bounds,
                  criterion.tol = 0,
                  initgrid.nx = 20, initgrid.nu = 20, audit.nx = 2500,
                  audit.nu = 25, audit.add = 100, audit.max = 25,
                  audit.tol,
                  point = FALSE, point.eyeweight = FALSE,
                  bootstraps = 0, bootstraps.m,
                  bootstraps.replace = TRUE,
                  levels = c(0.99, 0.95, 0.90), ci.type = 'backward',
                  specification.test = TRUE,
                  noisy = TRUE,
                  smallreturnlist = FALSE, seed = 12345, debug = FALSE) {
    call <- match.call(expand.dots = FALSE)
    envList <- list(m0 = environment(m0),
                    m1 = environment(m1),
                    ivlike = environment(ivlike),
                    parent = parent.frame())
    envProp <- try(environment(propensity), silent = TRUE)
    if (class(envProp) != "environment") {
        envList$propensity <- parent.frame()
    } else {
        envList$propensity <- envProp
    }
    for (i in 1:length(envList)) {
        if (is.null(envList[[i]])) envList[[i]] <- parent.frame()
    }

    ##---------------------------
    ## 1. Check linear programming dependencies
    ##---------------------------
    if (hasArg(lpsolver)) lpsolver <- tolower(lpsolver)
    if (is.null(lpsolver)) {
        if (requireNamespace("gurobi", quietly = TRUE)) {
            lpsolver <- "gurobi"
        } else if (requireNamespace("lpSolveAPI", quietly = TRUE)) {
            lpsolver <- "lpSolveAPI"
        } else if (requireNamespace("cplexAPI", quietly = TRUE)) {
            lpsolver <- "cplexAPI"
        } else {
            stop(gsub("\\s+", " ",
                      "Please install one of the following packages required for
                      estimation:
                      gurobi (version 7.5-1 or later);
                      cplexAPI (version 1.3.3 or later);
                      lpSolveAPI (version 5.5.2.0 or later)."),
                 call. = FALSE)
        }
    } else {
        if (! lpsolver %in% c("gurobi",
                              "cplexapi",
                              "lpsolveapi")) {
            stop(gsub("\\s+", " ",
                      paste0("Estimator is incompatible with linear programming
                             package '", lpsolver, "'. Please install one of the
                             following linear programming packages instead:
                             gurobi (version 7.5-1 or later);
                             cplexAPI (version 1.3.3 or later);
                             lpSolveAPI (version 5.5.2.0 or later).")),
                 call. = FALSE)
        }
    }
    if (hasArg(lpsolver.options)) {
        if (!is.list(lpsolver.options)) {
            stop(gsub("\\s+", " ",
                      paste0("'lpsolver.options' must be a list.
                               Each item in the list should correspond to an
                               option to be passed to the LP solver.
                               The name of the item should match the name
                               of the option, and the value of the item
                               should be the value to set the option to.")),
                 call. = FALSE)
        }
        if (hasArg(lpsolver.options.criterion) |
            hasArg(lpsolver.options.bounds)) {
            stop(gsub("\\s+", " ",
                      paste0("Either declare 'lpsolver.options'; or declare
                              'lpsolver.options.criterion' and/or
                              'lpsolver.options.bounds'; but not both.
                              In the case of
                              the latter, if only one set of options is
                              declared, then a set of default options will
                              be provided for the other.")),
                 call. = FALSE)
        }
    }
    if (hasArg(lpsolver.options.criterion)) {
        if (!is.list(lpsolver.options.criterion)) {
            stop(gsub("\\s+", " ",
                      paste0("'lpsolver.options.criterion' must be a list.
                               Each item in the list should correspond to an
                               option to be passed to the LP solver.
                               The name of the item should match the name
                               of the option, and the value of the item
                               should be the value to set the option to.")),
                 call. = FALSE)
        }
    }
    if (hasArg(lpsolver.options.bounds)) {
        if (!is.list(lpsolver.options.bounds)) {
            stop(gsub("\\s+", " ",
                      paste0("'lpsolver.options.bounds' must be a list.
                               Each item in the list should correspond to an
                               option to be passed to the LP solver.
                               The name of the item should match the name
                               of the option, and the value of the item
                               should be the value to set the option to.")),
                 call. = FALSE)
        }
    }
    if (hasArg(lpsolver.presolve)) {
        if (!is.logical(lpsolver.presolve)) {
            stop(paste0("'lpsolver.presolve' must either be TRUE or FALSE."),
                 call. = FALSE)
        }
        if (lpsolver != "gurobi") {
            warning(gsub("\\s+", " ",
                         paste0("The 'presolve' option is only implemented if
                                 the LP solver is Gurobi. For CPLEX and
                                 lp_solve, set the presolve parameter using
                                 'lpsolve.options', 'lpsolve.options.criterion',
                                 and 'lpsolve.options.bounds'.")),
                    call. = FALSE)
        }
        if ((hasArg(lpsolver.options) && !is.null(lpsolver.options$presolve)) |
            (hasArg(lpsolver.options.criterion) &&
             !is.null(lpsolver.options.criterion$presolve)) |
            (hasArg(lpsolver.options.bounds) &&
             !is.null(lpsolver.options.bounds$presolve))) {
            warning(gsub("\\s+", " ",
                         paste0("The 'presolve' option overrides the presolve
                                 parameters set in 'lpsolve.options',
                                 'lpsolve.options.criterion', and
                                 'lpsolve.options.bounds'.")),
                    call. = FALSE)
        }
    }
    if (debug) {
        if (lpsolver != "gurobi") {
            lpsolver <- "gurobi"
            warning(gsub("\\s+", " ",
                         "'debug = TRUE' is only permitted if
                          'lpsolver = \"gurobi\"'.
                           Linear programming output below is generated
                           by Gurobi."),
                    call. = FALSE, immediate. = TRUE)
        }
        if (!requireNamespace("gurobi", quietly = TRUE)) {
            stop(gsub("\\s+", " ",
                      "'debug = TRUE' is only permitted if
                       'lpsolver = \"gurobi\"'."))
        }
    }

    ##---------------------------
    ## 2. Check format of non-numeric arguments
    ##---------------------------
    ## Convert uname into a string
    uname <- deparse(substitute(uname))
    uname <- gsub("~", "", uname)
    uname <- gsub("\\\"", "", uname)
    ## Ensure data can be converted to a data.frame
    if (("data.frame" %in% class(data)) |
        ("matrix" %in% class(data))) {
        data <- as.data.frame(data)
    } else  {
        stop(gsub("\\s+", " ",
                  "'data' argument must either be a data.frame,
                   data.table, tibble, or matrix."),
             call. = FALSE)
    }
    ## Ensure MTRs are formulas
    if (classFormula(m0) && classFormula(m1)) {
        if (all(length(Formula::as.Formula(m0)) == c(0, 1)) &&
            all(length(Formula::as.Formula(m1)) == c(0, 1))) {
            mtrFail <- FALSE
        } else {
            mtrFail <- TRUE
        }
    } else {
        mtrFail <- TRUE
    }
    if (mtrFail) {
        stop(gsub("\\s+", " ",
                  "Arguments 'm0' and 'm1' must be one-sided formulas,
                   e.g. m0 = ~ 1 + u + x:I(u ^ 2)."),
             call. = FALSE)
    }
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
        tmpComp <- deparse(substitute(components))
        if (substr(tmpComp, 1, 2) != "l(" &&
            substr(tmpComp, 1, 2) == "c(") {
            stop(gsub("\\s+", " ",
                      "The 'components' argument should be declared
                       using 'l()' instead of 'c()'."),
                 call. = FALSE)
        }
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
                           data set and logical operators."),
                     call. = FALSE)
            } else if (subsetLogic == TRUE) {
                ## Currently prohobit logical vectors. Instead
                ## restrict to expressions only. Below is the code to
                ## revert to the case where logical vectors are
                ## allowed.
                ## subsetList <- list()
                ## subsetList[[1]] <- subset
                ## subset <- subsetList
                stop(gsub("\\s+", " ",
                          "Subset conditions should be logical
                           expressions involving variable names from the
                           data set and logical operators. Make sure the name of
                           the variable is not the same as that of another
                           object in the workspace. If so, rename the object
                           in the workspace, or the variable in the data set."),
                     call. = FALSE)
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
                           data set and logical operators."),
                         call. = FALSE)
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
                       conditions corresponds to using the full sample."),
                 call. = FALSE)
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
                      are logical.")),
                 call. = FALSE)
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
        treatStr <- deparse(substitute(treat))
        treatStr <- gsub("~", "", treatStr)
        treatStr <- gsub("\\\"", "", treatStr)
        if (! treatStr %in% colnames(data)) {
            stop("Declared treatment indicator not found in data",
                 call. = FALSE)
        }
    }
    if (hasArg(target)) {
        if (! target %in% c("ate", "att", "atu", "late", "genlate")) {
            stop(gsub("\\s+", " ",
                      "Specified target parameter is not recognized.
                      Choose from 'ate', 'att', 'atu', 'late', or 'genlate'."),
                 call. = FALSE)
        }
        if (target == "late") {
            ## Check the LATE arguments are declared properly.
            if (!(hasArg(late.to) & hasArg(late.from))) {
                stop(gsub("\\s+", " ",
                          "Target paramter 'late' requires arguments
                          'late.to', and 'late.from'."),
                     call. = FALSE)
            }
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
            late.Z = sort(names(late.from))
            late.Ztmp = sort(names(late.to))
            if (length(late.to) != length(late.from) |
                length(late.to) != length(late.Z)) {
                stop(gsub("\\s+", " ",
                          "The number of variables/values declared in 'late.to'
                           and 'late.from' must be equal."),
                     call. = FALSE)
            }
            if (!all(late.Z == late.Ztmp)) {
                stop(gsub("\\s+", " ",
                          "The variables declared in 'late.to' and 'late.from'
                           must be the same."),
                     call. = FALSE)
            }
            if (length(unique(late.Z)) != length(late.Z)) {
                stop(gsub("\\s+", " ",
                          "Each variable in 'late.to' and 'late.from' can only
                           be assigned one value."),
                     call. = FALSE)
            }
            if (!is.numeric(late.to) | !is.numeric(late.from)) {
                stop("Vectors 'late.from' and 'late.to' must be numeric.",
                     call. = FALSE)
            }
            late.from <- late.from[late.Z]
            late.to <- late.to[late.Z]
            if (all(late.to == late.from)) {
                stop(gsub("\\s+", " ",
                          "'late.to' must be different from 'late.from'."),
                     call. = FALSE)
            }
        }
        if (target == "genlate") {
            if (genlate.lb < 0 | genlate.ub > 1) {
                stop(gsub("\\s+", " ",
                          "'genlate.lb' and 'genlate.ub' must be between 0
                          and 1."), call. = FALSE)
            }
            if (genlate.lb >= genlate.ub) {
                stop("'genlate.lb' must be strictly less than 'genlate.ub'.",
                     call. = FALSE)
            }
        }
        if (target == "late" | target == "genlate") {
            ## Check that included insturments are declared properly
            if (hasArg(late.X)) {
                eval.X <- unlist(late.X)
                late.X <- names(late.X)
                if (!is.numeric(eval.X)) {
                    stop("Vector 'eval.X' must be numeric.",
                         call. = FALSE)
                }
            }
        }
        if (target != "genlate" &
            (hasArg("genlate.lb") | hasArg("genlate.ub"))) {
            warning(gsub("\\s+", " ",
                         "Unless target parameter is 'genlate', 'genlate.lb' and
                         'genlate.ub' arguments will not be used."))
        }
        if ((target != "late"  & target != "genlate") &
            (hasArg("eval.X") | hasArg("late.X"))) {
            warning(gsub("\\s+", " ",
                         "Unless target parameter is 'late' or
                          'genlate', the arguments eval.X' and
                          'late.X' will not be used."))
        }
        if ((hasArg(target.weight0) | hasArg(target.weight1))) {
            stop(gsub("\\s+", " ",
                      "A preset target weight is chosen, and a custom target
                      weight is provided. Please provide an input for 'target'
                      only, or for 'target.weight0' and 'target.weight1'."),
                 call. = FALSE)
        }
        target.weight0 <- NULL
        target.weight1 <- NULL
    } else {
        if (!(hasArg(target.weight0) & hasArg(target.weight1))) {
            stop(gsub("\\s+", " ",
                      "Only one target weight function is provided. If a custom
                      target weight is to be used, inputs for both
                      target.weight0 and target.weight1 must be provided."),
                 call. = FALSE)
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
                             ".")), call. = FALSE)
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
                             ".")), call. = FALSE)
        }
        if (length(target.weight0) != (length(target.knots0) + 1)) {
            stop(gsub("\\s+", " ",
                      paste0("The number of weight functions declared in
                                 target.weight0 must be exactly equal to the
                                 number of knots declared in target.knots0
                                 plus 1. Currently, the number of weights
                                 declared is ", length(target.weight0),
                             ", and the number of knots declared is ",
                             length(target.knots0), ".")),
                 call. = FALSE)
        }
        if (length(target.weight1) != (length(target.knots1) + 1)) {
            stop(gsub("\\s+", " ",
                      paste0("The number of weight functions declared in
                                 target.weight1 must be exactly equal to the
                                 number of knots declared in target.knots1
                                 plus 1. Currently, the number of weights
                                 declared is ", length(target.weight1),
                             ", and the number of knots declared is ",
                             length(target.knots1), ".")),
                 call. = FALSE)
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
                             ".")), call. = FALSE)
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
                             ".")), call. = FALSE)
            }
        }
    }
    if (hasArg(link)) {
        if (! link %in% c("linear", "logit", "probit")) {
            stop(gsub("\\s+", " ",
                      "Specified link is not recognized. Choose from 'linear',
                      'logit', or 'probit'."), call. = FALSE)
        }
    }
    if (hasArg(point)) {
        if (point == TRUE) {
            if (hasArg(m0.dec) | hasArg(m0.inc) |
                hasArg(m1.dec) | hasArg(m1.inc) |
                hasArg(mte.dec) | hasArg(mte.inc) |
                hasArg(audit.nu) | hasArg(audit.nx) | hasArg(audit.add) |
                hasArg(initgrid.nu) | hasArg(initgrid.nx)|
                hasArg(audit.max) | hasArg(audit.tol)) {
                warning(gsub("\\s+", " ",
                             "If argument 'point' is set to TRUE, then shape
                             restrictions on m0 and m1 are ignored, and the
                             audit procedure is not implemented."),
                        call. = FALSE)
            }
        }
    }
    ## Seed check
    if (!(is.numeric(seed) & length(seed) == 1)) {
        stop("'seed' must be a scalar.", call. = FALSE)
    }
    ## Audit checks
    if (!(is.numeric(criterion.tol) & criterion.tol >= 0)) {
        stop("Cannot set 'criterion.tol' below 0.", call. = FALSE)
    }
    if (!((initgrid.nu %% 1 == 0) & initgrid.nu >= 0)) {
        stop(gsub("\\s+", " ",
                  "'initgrid.nu' must be an integer than or equal to 0
                    (end points 0 and 1 are always included)."),
             call. = FALSE)
    }
    if (!((initgrid.nx %% 1 == 0) & initgrid.nx >= 0)) {
        stop("'initgrid.nx' must be an integer greater than or equal to 0.",
             call. = FALSE)
    }
    if (!((audit.nx %% 1 == 0) & audit.nx > 0)) {
        stop("'audit.nx' must be an integer greater than or equal to 1.",
             call. = FALSE)
    }
    if (audit.nx < initgrid.nx) {
        stop("'audit.nx' must be larger than 'initgrid.nx'.")
    }
    if (audit.nu < initgrid.nu) {
        stop("'audit.nu' must be larger than 'initgrid.nu'.")
    }
    if (!((audit.add %% 1 == 0) & audit.add > 0)) {
        stop("'audit.add' must be an integer greater than or equal to 1.",
             call. = FALSE)
    }
    if (hasArg(audit.tol)) {
        if (!is.numeric(audit.tol) |
            audit.tol < 0 |
            length(audit.tol) > 1) {
            stop("'audit.tol' must be a positive scalar.",
                 call. = FALSE)
        }
    }
    if (hasArg(m0.dec) | hasArg(m0.inc) |
        hasArg(m1.dec) | hasArg(m1.inc) |
        hasArg(mte.dec) | hasArg(mte.inc) |
        hasArg(m0.lb) | hasArg(m0.ub) |
        hasArg(m1.lb) | hasArg(m1.ub) |
        hasArg(mte.lb) | hasArg(mte.ub)) {
        noshape = FALSE ## indicator for whether shape restrictions declared
        if (!((audit.nu %% 1 == 0) & audit.nu >= 0)) {
            stop(gsub("\\s+", " ",
                      "'audit.nu' must be an integer than or equal to 0
                       (end points 0 and 1 are always included)."),
                 call. = FALSE)
        }
        if ((hasArg(m0.dec) && !is.logical(m0.dec)) |
            (hasArg(m1.dec) && !is.logical(m1.dec)) |
            (hasArg(mte.dec) && !is.logical(mte.dec)) |
            (hasArg(m0.inc) && !is.logical(m0.inc)) |
            (hasArg(m1.inc) && !is.logical(m1.inc)) |
            (hasArg(mte.inc) && !is.logical(mte.inc))) {
            stop(gsub("\\s+", " ",
                      "Monotonicity constraints 'm0.dec', 'm1.dec', 'mte.dec',
                       etc. must be either TRUE or FALSE."),
                 call. = FALSE)
        }
        if ((hasArg(m0.lb) && !is.numeric(m0.lb)) |
            (hasArg(m1.lb) && !is.numeric(m1.lb)) |
            (hasArg(mte.lb) && !is.numeric(mte.lb)) |
            (hasArg(m0.ub) && !is.numeric(m0.ub)) |
            (hasArg(m1.ub) && !is.numeric(m1.ub)) |
            (hasArg(mte.ub) && !is.numeric(mte.ub))) {
            stop(gsub("\\s+", " ",
                      "Boundedness constraints 'm0.lb', 'm1.lb', 'mte.lb',
                       etc. must be numeric."),
                 call. = FALSE)
        }
    } else {
        noshape = TRUE
        if (!((audit.nu %% 1 == 0) & audit.nu > 0)) {
            stop("'audit.nu' must be an integer greater than or equal to 1.",
                 call. = FALSE)
        }
    }
    if (!((audit.max %% 1 == 0) & audit.max > 0)) {
        stop("'audit.max' must be an integer greater than or equal to 1.",
             call. = FALSE)
    }
    ## Bootstrap checks
    if (bootstraps < 0 | bootstraps %% 1 != 0 | bootstraps == 1) {
        stop(gsub("\\s+", " ",
                  "'bootstraps' must either be 0, or be an integer greater than
                   or equal to 2."), call. = FALSE)
    }
    if (hasArg(bootstraps.m)) {
        if (bootstraps.m < 0 | bootstraps.m %% 1 != 0) {
            stop(gsub("\\s+", " ",
                      "'bootstraps.m' must be an integer greater than or equal
                      to 1."), call. = FALSE)
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
        stop("'bootstraps.replace' must be TRUE or FALSE.", call. = FALSE)
    }
    if (max(levels) >= 1 | min(levels) <= 0) {
        stop(gsub("\\s+", " ",
                  "'levels' must be a vector of values strictly between 0 and
                  1."), call. = FALSE)
    }
    levels <- sort(levels)
    if (! ci.type %in% c("forward", "backward", "both")) {
        stop(gsub("\\s+", " ",
                  "'ci.types' selects the type of confidence intervals to be
                  constructed for the treatment effect bound. It must be set to
                  either 'forward', 'backward', or 'both'."), call. = FALSE)
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
                      correctly."), call. = FALSE)
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
                          IV-like specifications."),
                     call. = FALSE)
            }
        }
    } else {
        stop(gsub("\\s+", " ",
                  "'ivlike' argument must either be a formula or a vector of
                  formulas."), call. = FALSE)
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
    parentFrame <- parent.frame()
    origm0 <- m0
    origm1 <- m1
    splinesobj <- list(removeSplines(m0, env = parentFrame),
                       removeSplines(m1, env = parentFrame))
    if (is.null(splinesobj[[1]]$formula)) {
        m0uCheck <- NULL
    } else {
        m0uCheck <- checkU(splinesobj[[1]]$formula, uname)
    }
    if (is.null(splinesobj[[2]]$formula)) {
        m1uCheck <- NULL
    } else {
        m1uCheck <- checkU(splinesobj[[2]]$formula, uname)
    }
    if (is.null(splinesobj[[1]]$splinesdict)) {
        m0uSplineCheck <- NULL
    } else {
        m0uSplineCheck <-
            checkU(as.formula(paste("~",
                                    paste(unlist(splinesobj[[1]]$splineslist),
                                          collapse = "+"))), uname)
    }
    if (is.null(splinesobj[[2]]$splinesdict)) {
        m1uSplineCheck <- NULL
    } else {
        m1uSplineCheck <-
            checkU(as.formula(paste("~",
                                    paste(unlist(splinesobj[[2]]$splineslist),
                                          collapse = "+"))), uname)
    }
    if (length(m0uCheck) + length(m0uSplineCheck) > 0) {
        message0 <- paste("  m0:",
                          paste(c(m0uCheck, m0uSplineCheck), collapse = ", "),
                          "\n")
    } else {
        message0 <- NULL
    }
    if (length(m1uCheck) + length(m1uSplineCheck) > 0) {
        message1 <- paste("  m1:",
                          paste(c(m1uCheck, m1uSplineCheck), collapse = ", "),
                          "\n")
    } else {
        message1 <- NULL
    }
    if (!is.null(message0) | !is.null(message1)) {
        if (uname != "x") {
            egv <- "x"
        } else {
            egv <- "v"
        }
        e1 <- paste0(uname, ", I(", uname, "^3)")
        e2 <- paste0(egv, ":", uname, ", ", egv, ":I(", uname, "^3)")
        e3 <- paste0("exp(", uname, "), I((", egv, " * ", uname, ")^2)")
        stop(gsub("\\s+", " ",
                  "The following terms are not declared properly."),
             "\n", message0, message1,
             gsub("\\s+", " ",
                  paste0("The unobserved variable '",
                         uname, "' must be declared as a
                         monomial, e.g. ", e1, ". The monomial can be
                         interacted with other variables, e.g. ",
                         e2, ". Expressions where the unobservable term
                         is not a monomial are either not permissable or
                         will not be parsed correctly,
                         e.g. ", e3, ". Try to rewrite the expression so
                         that '",  uname, "' is only included in monomials.")),
             call. = FALSE)
    }
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
                if (hasArg(target) && target == "late") {
                    if (!all(late.Z %in% vars_propensity)) {
                        stop (gsub("\\s+", " ",
                                   "All variables in 'late.to' and 'late.from'
                                    must be included in the propensity
                                    score model."), call. = FALSE)
                    }
                }
                if (hasArg(treat)) {
                    if (ptreat != treatStr) {
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
                          corresponds to propensity scores.")),
                         call. = FALSE)
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
                    treat <- treatStr
                }
                if (hasArg(target) && target == 'late') {
                    stop(gsub("\\s+", " ",
                              "Target parameter 'late' requires a propensity
                           score model. Agrument 'propensity' should be
                           a two-sided formula."),
                         call. = FALSE)
                }
            }
        } else {
            if (hasArg(target) && target == 'late') {
                stop(gsub("\\s+", " ",
                          "Target parameter 'late' requires a propensity
                           score model. Agrument 'propensity' should be
                           a two-sided formula."),
                     call. = FALSE)
            }
            propStr <- deparse(substitute(propensity))
            propStr <- gsub("~", "", propStr)
            propStr <- gsub("\\\"", "", propStr)
            if (!propStr %in% colnames(data)) {
                stop(gsub("\\s+", " ",
                          "Propensity score argument is interpreted as a
                          variable name, but is not found in the data set."),
                     call. = FALSE)
            }
            vars_propensity <- c(vars_propensity,
                                 propStr)
            ## Determine treatment variable
            if (hasArg(treat)) {
                treat <- treatStr
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
            propensity <- propStr
            propensity <- Formula::as.Formula(paste("~", propensity))
            if (length(all.vars(propensity)) > 1) {
                stop(gsub("\\s+", " ",
                          paste0("'propensity' argument must either be a
                          two-sided formula (if the propensity score is to be
                          estimated from the data), or a one-sided formula
                          containing a single variable on the RHS (where the
                          variable listed is included in the data, and
                          corresponds to propensity scores.")),
                     call. = FALSE)
            }
        }
    } else {
        ## Determine treatment variable
        if (hasArg(treat)) {
            treat <- treatStr
            vars_propensity <- treat
        } else if (is.list(ivlike)) {
            warning(gsub("\\s+", " ",
                         "First independent variable of first IV regression
                             is selected as the treatment variable."))
            treat <- all.vars(ivlike[[1]])[2]
        } else {
            stop("Treatment variable indeterminable.", call. = FALSE)
        }
        ## Construct propensity formula
        terms_propensity <- c(unlist(terms_formulas_x),
                              unlist(terms_formulas_z),
                              unlist(terms_mtr0),
                              unlist(terms_mtr1))
        ## Remove all u terms
        uterms <- c()
        um1 <- which(terms_propensity == uname)
        um2 <- grep(paste0("^[", uname, "][[:punct:]]"),
                    terms_propensity)
        um3 <- grep(paste0("[[:punct:]][", uname, "]$"),
                    terms_propensity)
        um4 <- grep(paste0("[[:punct:]][",
                           uname,
                           "][[:punct:]]"),
                    terms_propensity)
        um5 <- grep(paste0("[[:punct:]][", uname, "]\\s+"),
                    terms_propensity)
        um6 <- grep(paste0("\\s+[", uname, "][[:punct:]]"),
                    terms_propensity)
        um7 <- grep(paste0("\\s+[", uname, "]\\s+"),
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
                ## Remove all factor value specifications
                vars_components <- strsplit(vars_components, ":")
                vars_components <- lapply(vars_components,
                                          function(x) {
                                              xTmp <- sapply(x,
                                                             strsplit,
                                                             split = " - ")
                                              xTmp <- lapply(xTmp,
                                                             function(x) x[1])
                                              unlist(xTmp)
                                          })
                vars_components <- unique(unlist(vars_components))
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
    allvars <- allvars[allvars != uname]
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
                       in the propensity score model."),
                 call. = FALSE)
        }
        nLateX <- 0
        if (hasArg(late.X)) {
            nLateX <- length(late.X)
            if (!all(late.X %in% vars_propensity)) {
                stop(gsub("\\s+", " ",
                          "Variables in 'late.X' argument must be contained
                           in the propensity score model."),
                     call. = FALSE)
            }
        }
    }
    ## Keep only complete cases
    varFound1 <- sapply(allvars, exists, where = data)
    vars_data <- allvars[varFound1]
    missingVars <- allvars[!varFound1]
    missingVars <- missingVars[!missingVars %in% c("intercept", uname)]
    for (i in 1:length(envList)) {
        varFound2 <- sapply(missingVars, exists, where = envList[[i]])
        missingVars <- missingVars[!varFound2]
    }
    if (length(missingVars) > 0) {
        varError <- paste0("The following variables are not contained
                          in the data set or in the parent.frame(): ",
                          paste(missingVars, collapse = ", "),
                          ".")
        stop(gsub("\\s+", " ", varError), call. = FALSE)
    }
    data  <- data[complete.cases(data[, vars_data]), ]
    ## Adjust row names to handle bootstrapping
    rownames(data) <- as.character(seq(1, nrow(data)))
    ## Construct a list of what variables interact with the
    ## spline.
    splinesobj <- interactSplines(splinesobj, origm0, origm1, data, uname)
    ## Check that all boolean variables have non-zero variance, and
    ## that all factor variables are complete.
    allterms <- unlist(c(terms_formulas_x, terms_formulas_z, terms_mtr0,
                         terms_mtr1, terms_components))
    allterms <- unique(unlist(sapply(allterms, strsplit, split = ":")))
    allterms <- parenthBoolean(allterms)
    ## Check factors
    factorPos <- grep("factor\\([[:alnum:]]*\\)", allterms)
    if (length(factorPos) > 0) {
        factorDict <- list()
        for (i in allterms[factorPos]) {
            var <- substr(i, start = 8, stop = (nchar(i) - 1))
            factorDict[[var]] <- sort(unique(data[, var]))
        }
    } else {
        factorDict <- NULL
    }
    ## To check for binary variables
    binMinCheck <- apply(data, 2, min)
    binMaxCheck <- apply(data, 2, max)
    binUniqueCheck <- apply(data, 2, function(x) length(unique(x)))
    binaryPos <- as.logical((binMinCheck == 0) *
                           (binMaxCheck == 1) *
                           (binUniqueCheck == 2))
    binaryVars <- colnames(data)[binaryPos]
    if (length(binaryVars) == 0) binaryVars <- NULL
    ## To check booleans expressions, you just need to ensure there is
    ## non-zero variance in the design matrix.
    boolPos <- NULL
    for (op in c("==", "!=", ">=", "<=", ">",  "<")) {
        boolPos <- c(boolPos, grep(op, allterms))
    }
    if (length(boolPos) > 0) {
        boolPos <- unique(boolPos)
        boolVars <- allterms[boolPos]
    } else {
        boolVars <- NULL
    }
    ## Save call options for return
    opList <- modcall(call,
                      newcall = list,
                      keepargs =  c('target', 'late.from', 'late.to',
                                    'late.X', 'genlate.lb',
                                    'genlate.ub', 'm0', 'm1', 'm1.ub',
                                    'm0.ub', 'm1.lb', 'm0.lb',
                                    'mte.ub', 'mte.lb', 'm0.dec',
                                    'm0.inc', 'm1.dec', 'm1.inc',
                                    'mte.dec', 'mte.inc', 'link',
                                    'seed'))
    opList <- eval(opList)
    opList$ivlike <- ivlike
    opList$components <- components
    if (hasArg(subset)) opList$subset <- subset
    if (hasArg(propensity)) opList$propensity <- propensity
    if (hasArg(uname)) opList$uname <- uname
    if (hasArg(treat)) opList$treat <- treat
    if (hasArg(target.weight0)) opList$target.weight0 <- target.weight0
    if (hasArg(target.weight1)) opList$target.weight1 <- target.weight1
    if (hasArg(target.knots0)) opList$target.knots0 <- target.knots0
    if (hasArg(target.knots1)) opList$target.knots1 <- target.knots1

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
                                         "late.X", "eval.X",
                                         "specification.test"),
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
                                           vars_data = quote(vars_data),
                                           terms_mtr0 = quote(terms_mtr0),
                                           terms_mtr1 = quote(terms_mtr1),
                                           treat = quote(treat),
                                           propensity = quote(propensity),
                                           splinesobj = quote(splinesobj),
                                           components = quote(components),
                                           environments = quote(envList)))
    if (hasArg(target) && (target == "late" | target == "genlate")) {
        if (target == "late") {
            estimateCall <- modcall(estimateCall,
                                    newargs = list(late.Z = late.Z,
                                                   late.to = late.to,
                                                   late.from = late.from))
        }
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
            if (!noisy) {
                ## Some output must be returned, evne if noisy = FALSE
                cat("\n")
                cat("Bounds on the target parameter: [",
                    fmtResult(origEstimate$bounds[1]), ", ",
                fmtResult(origEstimate$bounds[2]), "]\n",
                sep = "")
                if (origEstimate$audit.count == 1) rs <- "round.\n"
                if (origEstimate$audit.count > 1) rs <- "rounds.\n"
                if (origEstimate$audit.count < audit.max) {
                    cat("Audit terminated successfully after",
                        origEstimate$audit.count,
                        rs, "\n")
                }
                if (origEstimate$audit.count == audit.max) {
                    if (is.null(origEstimate$audit.grid$violations)) {
                        cat("Audit terminated successfully after",
                            origEstimate$audit.count,
                            rs, "\n")
                    } else {
                        cat("\n")
                        warning(gsub("\\s+", " ",
                                     "Audit reached audit.max.
                                      Try increasing audit.max."),
                                call. = FALSE,
                                immediate. = TRUE)
                    }
                }
            }
            output <- origEstimate
            output$call.options <- opList
            output <- output[sort(names(output))]
            class(output) <- "ivmte"
            return(invisible(output))
        } else {
            ## Obtain audit grid from original estimate
            audit.grid <- list(support = origEstimate$audit.grid$audit.x,
                               uvec = origEstimate$audit.grid$audit.u)
            if (is.null(audit.grid$support)) audit.grid$noX <- TRUE
            if (!is.null(audit.grid$support)) audit.grid$noX <- FALSE
            ## Estimate bounds with resampling
            set.seed(seed)
            bseeds <- round(runif(bootstraps) * 1000000)
            boundEstimates <- NULL
            propEstimates <- NULL
            b <- 1
            bootCriterion <- NULL
            bootFailN <- 0
            bootFailNote <- ""
            bootFailIndex <- NULL
            ## Turn off specification test if criterion is 0
            if (origEstimate$audit.criterion == 0) {
                specification.test <- FALSE
            }
            if (specification.test) {
                origSset <- lapply(origEstimate$sset, function(x) {
                    x[c("ivspec", "beta", "g0", "g1")]
                })
                origCriterion <- origEstimate$audit.criterion
            } else {
                origSset <- NULL
                origCriterion <- NULL
            }
            if (!hasArg(bootstraps.m)) bootstraps.m <- nrow(data)
            if (bootstraps.m > nrow(data) && bootstraps.replace == FALSE) {
                stop(gsub("\\s+", " ",
                          "Cannot draw more observations than the number of rows
                           in the data set when 'bootstraps.replace = FALSE'."),
                     call. = FALSE)
            }
            while (b <= bootstraps) {
                set.seed(bseeds[b])
                bootIDs  <- sample(x = seq(1, nrow(data)),
                                   size = bootstraps.m,
                                   replace = bootstraps.replace)
                bdata <- data[bootIDs, ]
                if (noisy == TRUE) {
                    cat("Bootstrap iteration ", b, "...\n", sep = "")
                }
                bootCall <-
                    modcall(estimateCall,
                            dropargs = c("data", "noisy", "seed",
                                         "audit.grid",
                                         "count.moments"),
                            newargs = list(data = quote(bdata),
                                           noisy = FALSE,
                                           seed = bseeds[b],
                                           audit.grid = audit.grid,
                                           orig.sset = origSset,
                                           orig.criterion = origCriterion,
                                           count.moments = FALSE))
                bootEstimate <- try(eval(bootCall), silent = TRUE)
                if (is.list(bootEstimate)) {
                    boundEstimates  <- rbind(boundEstimates,
                                             bootEstimate$bound)
                    if (specification.test) {
                        bootCriterion <- c(bootCriterion,
                                           bootEstimate$specification.test)
                    }
                    if (!"propensity.coef" %in% names(bootEstimate) &&
                        (class(bootEstimate$propensity$model)[1] == "lm" |
                         class(bootEstimate$propensity$model)[1] == "glm")) {
                        propEstimates <-
                            cbind(propEstimates,
                                  bootEstimate$propensity$model$coef)
                    } else {
                        propEstimates <- cbind(propEstimates,
                                               bootEstimate$propensity.coef)
                    }
                    b <- b + 1
                    if (noisy == TRUE) {
                        cat("    Audit count: ",
                            bootEstimate$audit.count, "\n", sep = "")
                        cat("    Minimum criterion: ",
                            fmtResult(bootEstimate$audit.criterion),
                            "\n", sep = "")
                        cat("    Bounds: ",
                            paste0("[",
                                   fmtResult(bootEstimate$bounds[1]),
                                   ", ",
                                   fmtResult(bootEstimate$bounds[2]),
                                   "]"), "\n\n", sep = "")
                    }
                } else {
                    if (noisy) {
                        warning(paste0(bootEstimate, ", resampling...\n"),
                                call. = FALSE,
                                immediate. = TRUE)
                    }
                    bseeds[b] <- round(runif(1) * 1000000)
                    bootFailN <- bootFailN + 1
                    bootFailIndex <- unique(c(bootFailIndex, b))
                }
            }
            if (noisy) {
                cat("--------------------------------------------------\n")
                cat("Results", "\n")
                cat("--------------------------------------------------\n")
            }
            cat("\n")
            ## Some output must be returned, evne if noisy = FALSE
            cat("Bounds on the target parameter: [",
                fmtResult(origEstimate$bounds[1]), ", ",
                fmtResult(origEstimate$bounds[2]), "]\n",
                sep = "")
            if (origEstimate$audit.count == 1) rs <- "round."
            if (origEstimate$audit.count > 1) rs <- "rounds."
            if (origEstimate$audit.count < audit.max) {
                cat("Audit terminated successfully after",
                    origEstimate$audit.count,
                    rs, "\n")
            }
            if (origEstimate$audit.count == audit.max) {
                if (is.null(origEstimate$audit.grid$violations)) {
                    cat("Audit terminated successfully after",
                        origEstimate$audit.count,
                        rs, "\n")
                } else {
                    warning(gsub("\\s+", " ",
                                 "Audit reached audit.max.
                                      Try increasing audit.max."),
                            call. = FALSE,
                            immediate. = TRUE)
                }
            }
            ## Obtain standard errors of bounds
            bootSE <- apply(boundEstimates, 2, sd)
            ## Construct confidence intervals for bounds
            ci <- boundCI(bounds = origEstimate$bounds,
                          bounds.resamples = boundEstimates,
                          n = nrow(data),
                          m = bootstraps.m,
                          levels = levels,
                          type = "both")
            ## Construct confidence intervals for propensity scores
            if (!is.null(propEstimates)) {
                propensity.ci <- list()
                if (!"propensity.coef" %in% names(origEstimate) &&
                    class(origEstimate$propensity$model)[1]
                    != "data.frame" &&
                       class(origEstimate$propensity$model)[1]
                    != "data.table") {
                    propCoef <- origEstimate$propensity$model$coef
                } else {
                    propCoef <- origEstimate$propensity.coef
                }
                propSE  <- apply(propEstimates, 1, sd)
                for (level in levels) {
                    pLower <- (1 - level) / 2
                    pUpper <- 1 - (1 - level) / 2
                    probVec <- c(pLower, pUpper)
                    tmpPropCi1 <- apply(propEstimates, 1, quantile,
                                        probs = probVec,
                                        type = 1)
                    tmpPropCi2 <- sweep(x = tcrossprod(c(qnorm(pLower),
                                                         qnorm(pUpper)),
                                                       propSE), MARGIN = 2,
                                        propCoef, FUN = "+")
                    colnames(tmpPropCi2) <- colnames(tmpPropCi1)
                    rownames(tmpPropCi2) <- rownames(tmpPropCi1)
                    propensity.ci$ci1[[paste0("level", level * 100)]] <-
                        t(tmpPropCi1)
                    propensity.ci$ci2[[paste0("level", level * 100)]] <-
                        t(tmpPropCi2)
                }
            }
            ## Obtain p-value
            pvalue <- c(boundPvalue(bounds = origEstimate$bounds,
                                    bounds.resamples = boundEstimates,
                                    n = nrow(data),
                                    m = bootstraps.m,
                                    type = "backward"),
                        boundPvalue(bounds = origEstimate$bounds,
                                    bounds.resamples = boundEstimates,
                                    n = nrow(data),
                                    m = bootstraps.m,
                                    type = "forward"))
            names(pvalue) <- c("backward", "forward")
            if (ci.type == "both")  {
                ciTypes <- c("backward", "forward")
            } else {
                ciTypes <- ci.type
            }
            ciN <- 1
            for (i in ciTypes) {
                cat("\nBootstrapped confidence intervals (",
                    i, "):\n", sep = "")
                for (j in 1:length(levels)) {
                    cistr <- paste0("[",
                                    fmtResult(ci[[i]][j, 1]),
                                    ", ",
                                    fmtResult(ci[[i]][j, 2]),
                                    "]")
                    cat("    ",
                        levels[j] * 100,
                        "%: ",
                        cistr, "\n", sep = "")
                }
                cat("p-value: ", fmtResult(pvalue[ciN]), "\n", sep = "")
                ciN <- ciN + 1
            }
            ## Obtain specification test
            if (specification.test) {
                criterionPValue <- mean(origEstimate$audit.criterion <=
                                        bootCriterion)
                cat("Bootstrapped specification test p-value: ",
                    criterionPValue, "\n", sep = "")
            }
            ## Print bootstrap statistics
            cat(sprintf("Number of bootstraps: %s",
                        bootstraps))
            if (bootFailN > 0) {
                cat(sprintf(" (%s failed and redrawn)\n",
                            bootFailN))
            } else {
                cat("\n")
            }
            cat("\n")
            ## Return output
            output <- c(origEstimate,
                        list(bounds.se = bootSE,
                             bounds.bootstraps = boundEstimates,
                             bounds.ci = ci,
                             pvalue = pvalue,
                             bootstraps = bootstraps,
                             bootstraps.failed = bootFailN))
            if (specification.test) {
                output$specification.pvalue <- criterionPValue
            }
            if (!is.null(propEstimates)) {
                output$propensity.se <- propSE
                output$propensity.ci  <- propensity.ci
            }
            output$call.options <- opList
            output <- output[sort(names(output))]
            class(output) <- "ivmte"
            return(invisible(output))
        }
    }
    ## Point estimate without resampling
    if (point == TRUE & bootstraps == 0) {
        origEstimate <- eval(estimateCall)
        if (!noisy) {
            ## Some output must be returned, even if noisy = FALSE
            cat("\n")
            cat("Point estimate of the target parameter: ",
                fmtResult(origEstimate$pointestimate), "\n\n",
                sep = "")
        }
        output <- origEstimate
        output$call.options = opList
        output <- output[sort(names(output))]
        class(output) <- "ivmte"
        return(invisible(output))
    }
    ## Point estimate with resampling
    if (point == TRUE & bootstraps > 0) {
        set.seed(seed)
        bseeds <- round(runif(bootstraps) * 1000000)
        origEstimate <- eval(estimateCall)
        teEstimates  <- NULL
        mtrEstimates <- NULL
        propEstimates <- NULL
        jstats <- NULL
        b <- 1
        bootFailN <- 0
        bootFailNote <- ""
        bootFailIndex <- NULL
        if (!hasArg(bootstraps.m)) bootstraps.m <- nrow(data)
        if (bootstraps.m > nrow(data) && bootstraps.replace == FALSE) {
            stop(gsub("\\s+", " ",
                      "Cannot draw more observations than the number of rows
                           in the data set when 'bootstraps.replace = FALSE'."),
                 call. = FALSE)
        }
        totalBootstraps <- 0
        factorMessage <- 0
        factorText <- gsub("\\s+", " ",
                           "Insufficient variation in categorical variables
                           (i.e. factor variables, binary variables,
                           boolean expressions) in the bootstrap sample.
                           Additional bootstrap samples will be drawn.")
        while (b <= bootstraps) {
            set.seed(bseeds[b])
            totalBootstraps <- totalBootstraps + 1
            bootIDs  <- sample(seq(1, nrow(data)),
                                 size = bootstraps.m,
                               replace = bootstraps.replace)
            bdata <- data[bootIDs, ]
            ## Check if the bootstrap data contains sufficient
            ## variation in all boolean and factor expressions.
            if (!is.null(factorDict)) {
                for (i in 1:length(factorDict)) {
                    factorCheck <-
                        suppressWarnings(
                            all(sort(unique(bdata[, names(factorDict)[i]])) ==
                            factorDict[[i]]))
                    if (!factorCheck) break
                }
                if (!factorCheck) {
                    if (factorMessage == 0) {
                        factorMessage <- 1
                        warning(factorText, call. = FALSE, immediate. = TRUE)
                    }
                    next
                    }
            }
            ## Check if the binary variables contain sufficient variation
            if (!is.null(binaryVars)) {
                binarySd <- apply(as.matrix(bdata[, binaryVars]), 2, sd)
                if (! all(binarySd > 0)) {
                    if (factorMessage == 0) {
                        factorMessage <- 1
                        warning(factorText, call. = FALSE, immediate. = TRUE)
                    }
                    next
                    }
            }
            ## Check if the boolean variables contain sufficient variation
            if (!is.null(boolVars)) {
                boolFormula <-
                    as.formula(paste("~ 0 +", paste(boolVars,
                                                    collapse = " + ")))
                boolDmat <- design(boolFormula, bdata)$X
                boolSd <- apply(boolDmat, 2, sd)
                if (! all(boolSd > 0)) {
                    if (factorMessage == 0) {
                        factorMessage <- 1
                        warning(factorText, call. = FALSE, immediate. = TRUE)
                    }
                    next
                }
            }
            if (noisy == TRUE) {
                cat("Bootstrap iteration ", b, "...\n", sep = "")
            }
            bootCall <-
                modcall(estimateCall,
                        dropargs = c("data", "noisy", "seed"),
                        newargs = list(data = quote(bdata),
                                       noisy = FALSE,
                                       seed = bseeds[b],
                                       point.center = origEstimate$moments,
                                       point.redundant =
                                           origEstimate$redundant))
            bootEstimate <- try(eval(bootCall), silent = TRUE)
            if (is.list(bootEstimate)) {
                teEstimates  <- c(teEstimates, bootEstimate$pointestimate)
                mtrEstimates <- cbind(mtrEstimates, bootEstimate$mtr.coef)
                if (!"propensity.coef" %in% names(bootEstimate)) {
                    propEstimates <- cbind(propEstimates,
                                           bootEstimate$propensity$model$coef)
                } else {
                    propEstimates <- cbind(propEstimates,
                                           bootEstimate$propensity.coef)
                }
                if (!is.null(bootEstimate$jtest)) {
                    jstats <- c(jstats, bootEstimate$jtest[1])
                }
                b <- b + 1
                if (noisy == TRUE) {
                    cat("    Point estimate:",
                        fmtResult(bootEstimate$pointestimate), "\n\n", sep = "")
                }
            } else {
                if (noisy == TRUE) {
                    message("    Error, resampling...\n", sep = "")
                }
                bseeds[b] <- round(runif(1) * 1000000)
                bootFailN <- bootFailN + 1
                bootFailIndex <- unique(c(bootFailIndex, b))
            }
        }
        if (bootFailN > 0) {
            warning(gsub("\\s+", " ",
                         paste0("Bootstrap iteration(s) ",
                                paste(bootFailIndex, collapse = ", "),
                                " failed. Failed bootstraps are repeated.")),
                    call. = FALSE)
        }
        bootSE <- sd(teEstimates)
        mtrSE  <- apply(mtrEstimates, 1, sd)
        if (!is.null(propEstimates)) {
            propSE  <- apply(propEstimates, 1, sd)
        }
        ## Construct p-values (point estimate and J-test)
        pvalue <- c(nonparametric =
                        (sum(teEstimates - origEstimate$pointestimate >=
                             abs(origEstimate$pointestimate)) +
                         sum(teEstimates - origEstimate$pointestimate <=
                             -abs(origEstimate$pointestimate))) / bootstraps,
                    parametric =
                        pnorm(-abs(origEstimate$pointestimate -
                                   mean(teEstimates)) / bootSE) * 2)
        if (!is.null(jstats)) {
            jtest <- c(mean(jstats >= origEstimate$jtest[1]),
                       origEstimate$jtest)
            names(jtest) <- c("Bootstrapped p-value", names(origEstimate$jtest))
            jtest <- jtest[c(2, 1, 3, 4)]
        } else {
            jtest <- NULL
        }
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
            if (!is.null(propEstimates)) {
                assign(paste0("propci1", level * 100),
                       apply(propEstimates, 1, quantile,
                             probs = probVec,
                             type = 1))
            }
            ## Conf. int. 2: percentile method using Z statistics
            tmpCi2 <- origEstimate$pointestimate +
                c(qnorm(pLower), qnorm(pUpper)) * bootSE
            names(tmpCi2) <- paste0(probVec * 100, "%")
            tmpMtrCi2 <- sweep(x = tcrossprod(c(qnorm(pLower),
                                                qnorm(pUpper)),
                                              mtrSE),
                               MARGIN = 2, origEstimate$mtr.coef, FUN = "+")
            colnames(tmpMtrCi2) <- colnames(get(paste0("mtrci1",
                                                       level * 100)))
            rownames(tmpMtrCi2) <- rownames(get(paste0("mtrci1",
                                                       level * 100)))
            assign(paste0("ci2", level * 100), tmpCi2)
            assign(paste0("mtrci2", level * 100), tmpMtrCi2)
            ## Special case for propensity scores, which may not exist
            if (!"propensity.coef" %in% names(origEstimate)) {
                propCoef <- origEstimate$propensity$model$coef
            } else {
                propCoef <- origEstimate$propensity.coef
            }
            if (!is.null(propEstimates)) {
                tmpPropCi2 <- sweep(x = tcrossprod(c(qnorm(pLower),
                                                     qnorm(pUpper)),
                                                   propSE), MARGIN = 2,
                                    propCoef, FUN = "+")
                colnames(tmpPropCi2) <- colnames(get(paste0("propci1",
                                                            level * 100)))
                rownames(tmpPropCi2) <- rownames(get(paste0("propci1",
                                                            level * 100)))
                assign(paste0("propci2", level * 100), tmpPropCi2)
            }
        }
        ## Prepare output
        output1 <- c(origEstimate,
                     list(pointestimate.se = bootSE,
                          mtr.se = mtrSE))
        if (!is.null(propEstimates)) output1$propensity.se <- propSE
        output1 <- c(output1,
                     list(pointestimate.bootstraps = teEstimates,
                          mtr.bootstraps = t(mtrEstimates)))
        pointestimate.ci <- list()
        mtr.ci <- list()
        if (!is.null(propEstimates)) {
            propensity.ci <- list()
        }
        for (level in levels) {
            pointestimate.ci$nonparametric <-
                rbind(pointestimate.ci$nonparametric,
                      get(paste0("ci1", level * 100)))
            pointestimate.ci$normal <-
                rbind(pointestimate.ci$normal,
                      get(paste0("ci2", level * 100)))
            mtr.ci$nonparametric[[paste0("level", level * 100)]] <-
                t(get(paste0("mtrci1", level * 100)))
            mtr.ci$normal[[paste0("level", level * 100)]] <-
                t(get(paste0("mtrci2", level * 100)))
            if (!is.null(propEstimates)) {
                propensity.ci$nonparametric[[paste0("level", level * 100)]] <-
                    t(get(paste0("propci1", level * 100)))
                propensity.ci$normal[[paste0("level", level * 100)]] <-
                    t(get(paste0("propci2", level * 100)))
            }
        }
        rownames(pointestimate.ci$nonparametric) <- levels
        rownames(pointestimate.ci$normal) <- levels
        colnames(pointestimate.ci$nonparametric) <- c("lb", "ub")
        colnames(pointestimate.ci$normal) <- c("lb", "ub")
        output2 <- list(pointestimate.ci = pointestimate.ci,
                        mtr.ci = mtr.ci)
        if (!is.null(propEstimates)) {
            output2$propensity.ci = propensity.ci
        }
        output3 <- list(pvalue = pvalue,
                        bootstraps = bootstraps,
                        bootstraps.failed = bootFailN,
                        jtest = jtest,
                        jtest.bootstraps = jstats)
        if ("jtest" %in% names(output1) &&
            "jtest" %in% names(output3)) {
            output1$jtest <- NULL
        }
        output <- c(output1, output2, output3)
        if (noisy) {
            cat("--------------------------------------------------\n")
            cat("Results", "\n")
            cat("--------------------------------------------------\n")
        }
        cat("\nPoint estimate of the target parameter: ",
            fmtResult(origEstimate$pointestimate), "\n",
            sep = "")
        cat("\nBootstrapped confidence intervals (nonparametric):\n")
        for (level in levels) {
            ci1str <- get(paste0("ci1", level * 100))
            ci1str <- paste0("[",
                             fmtResult(ci1str[1]),
                             ", ",
                             fmtResult(ci1str[2]),
                             "]")
            cat("    ",
                level * 100,
                "%: ",
                ci1str, "\n", sep = "")
        }
        cat("p-value: ",
            fmtResult(pvalue[1]), "\n", sep = "")
        if (!is.null(jtest)) {
            cat("Bootstrapped J-test p-value: ",
                fmtResult(jtest[2]), "\n", sep = "")
        }
        cat(sprintf("Number of bootstraps: %s",
                    bootstraps))
        if (bootFailN > 0) {
            cat(sprintf(" (%s failed and redrawn)\n",
                        bootFailN))
        } else {
            cat("\n")
        }
        cat("\n")
        ## cat("Bootstrapped confidence intervals (normal quantiles):\n")
        ## for (level in levels) {
        ##     ci2str <- get(paste0("ci2", level * 100))
        ##     ci2str <- paste0("[",
        ##                      fmtResult(ci2str[1]),
        ##                      ", ",
        ##                      fmtResult(ci2str[2]),
        ##                      "]")
        ##     cat("    ",
        ##         level * 100,
        ##         "%: ",
        ##         ci2str, "\n", sep = "")
        ## }
        ## cat("p-value: ",
        ##     fmtResult(pvalue[2]), "\n\n", sep = "")
        if (totalBootstraps > bootstraps) {
            warning(gsub("\\s+", " ",
                         paste0("In order to obtain ", bootstraps, " boostrap
                         samples without omiting any
                         levels from all categorical variables,
                         a total of ", totalBootstraps, " samples
                         had to be drawn. This is due to factor variables
                         and/or boolean variables potentially being omitted from
                         boostrap samples.")), call. = FALSE)
        }
        output$call.options <- opList
        output <- output[sort(names(output))]
        class(output) <- "ivmte"
        return(invisible(output))
    }
}

#' Construct confidence intervals for treatment effects under partial
#' identification
#'
#' This function constructs the forward and backward confidence
#' intervals for the treatment effect under partial identification.
#'
#' @param bounds vector, bounds of the treatment effects under partial
#'     identification.
#' @param bounds.resamples matrix, stacked bounds of the treatment
#'     effects under partial identification. Each row corresponds to a
#'     subset resampled from the original data set.
#' @param n integer, size of original data set.
#' @param m integer, size of resampled data sets.
#' @param levels vector, real numbers between 0 and 1. Values
#'     correspond to the level of the confidence intervals constructed
#'     via bootstrap.
#' @param type character. Set to 'forward' to construct the forward
#'     confidence interval for the treatment effect bounds. Set to
#'     'backward' to construct the backward confidence interval for
#'     the treatment effect bounds. Set to 'both' to construct both
#'     types of confidence intervals.
#' @return if \code{type} is 'forward' or 'backward', then the
#'     corresponding type of confidence interval for each level is
#'     returned. The output is in the form of a matrix, with each row
#'     corresponding to a level. If \code{type} is 'both', then a list
#'     is returned. One element of the list is the matrix of backward
#'     confidence intervals, and the other element of the list is the
#'     matrix of forward confidence intervals.
boundCI <- function(bounds, bounds.resamples, n, m, levels, type) {
    if (type == "both") output <- list()
    ## Rescale and center bounds
    boundLBmod <- sqrt(m) *
        (bounds.resamples[, 1] - bounds[1])
    boundUBmod <- sqrt(m) *
        (bounds.resamples[, 2] - bounds[2])
    ## Construct backward confidence interval
    if (type %in% c('backward', 'both')) {
        backwardCI <- cbind(bounds[1] +
                            (1 / sqrt(n)) *
                            quantile(boundLBmod,
                                     0.5 * (1 - levels),
                                     type = 1),
                            bounds[2] +
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
        forwardCI <- cbind(bounds[1] -
                           (1 / sqrt(n)) *
                           quantile(boundLBmod,
                                    0.5 * (1 + levels),
                                    type = 1),
                           bounds[2] -
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
#' @param bounds vector, bounds of the treatment effects under partial
#'     identification.
#' @param bounds.resamples matrix, stacked bounds of the treatment
#'     effects under partial identification. Each row corresponds to a
#'     subset resampled from the original data set.
#' @param n integer, size of original data set.
#' @param m integer, size of resampled data sets.
#' @param type character. Set to 'forward' to construct the forward
#'     confidence interval for the treatment effect bounds. Set to
#'     'backward' to construct the backward confidence interval for
#'     the treatment effect bounds. Set to 'both' to construct both
#'     types of confidence intervals.
#' @return If \code{type} is 'forward' or 'backward', a scalar p-value
#'     corresponding to the type of confidence interval is
#'     returned. If \code{type} is 'both', a vector of p-values
#'     corresponding to the forward and backward confidence intervals
#'     is returned.
boundPvalue <- function(bounds, bounds.resamples, n, m, type) {
    levels <- seq(0, 1, 1 / nrow(bounds.resamples))[-1]
    cis <- boundCI(bounds, bounds.resamples, n, m, levels, type)
    posNo0 <- apply(cis, MARGIN = 1, FUN = function(x) {
        !((0 >= x[1]) && (0 <= x[2]))
    })
    if (all(posNo0)) {
        return(0)
    } else if (all(!posNo0)) {
        return(1)
    } else {
        return(1 - max(levels[posNo0]))
    }
}

#' Check polynomial form of the u-term
#'
#' This function ensures that the unobservable term enters into the
#' MTR in the correct manner. That is, it enters as a polynomial.
#'
#' @param formula a formula.
#' @param uname name of the unobserved variable.
#' @return If the unobservable term is entered correctly into the
#'     formula, then \code{NULL} is returned. Otherwise, the vector of
#'     incorrect terms is returned.
checkU <- function(formula, uname) {
    termsList <- attr(terms(formula), "term.labels")
    termsList <- unique(unlist(strsplit(termsList, ":")))
    termsVarList <- lapply(termsList, function(x) {
        all.vars(as.formula(paste("~", x)))
    })
    upos <- unlist(lapply(termsVarList, function(x) uname %in% x))
    parPos <- unlist(lapply(termsList, function(x) grepl("\\(", x)))
    termsList <- termsList[as.logical(upos * parPos)]
    checkVec <- grepl(paste0("^I\\(", uname, "\\^[0-9]+\\)$"), termsList)
    if (all(checkVec))  errorTermsFormula <- NULL
    if (!all(checkVec)) errorTermsFormula <- termsList[!checkVec]
    return(errorTermsFormula)
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
#' @param late.Z vector of variable names used to define the LATE.
#' @param late.from baseline set of values of Z used to define the
#'     LATE.
#' @param late.to comparison set of values of Z used to define the
#'     LATE.
#' @param late.X vector of variable names of covariates which we
#'     condition on when defining the LATE.
#' @param eval.X numeric vector of the values at which we condition
#'     variables in \code{late.X} on when estimating the LATE.
#' @param audit.grid list, contains the A A matrix used in the audit
#'     for the original sample, as well as the RHS vector used in the
#'     audit from the original sample.
#' @param point.center numeric, a vector of GMM moment conditoins
#'     evaluated at a solution. When bootstrapping, the moment
#'     conditions from the original sample can be passed through this
#'     argument to recenter the bootstrap distribution of the
#'     J-statistic.
#' @param point.redundant vector of integers indicating which
#'     components in the S-set are redundant.
#' @param count.moments boolean, indicate if number of linearly
#'     independent moments should be counted.
#' @param orig.sset list, only used for bootstraps. The list caontains
#'     the gamma moments for each element in the S-set, as well as the
#'     IV-like coefficients.
#' @param orig.criterion numeric, only used for bootstraps. The scalar
#'     corresponds to the minimum observational equivalence criterion
#'     from the original sample.
#' @param vars_y character, variable name of observed outcome
#'     variable.
#' @param vars_mtr character, vector of variables entering into
#'     \code{m0} and \code{m1}.
#' @param terms_mtr0 character, vector of terms entering into
#'     \code{m0}.
#' @param terms_mtr1 character, vector of terms entering into
#'     \code{m1}.
#' @param vars_data character, vector of variables that can be found
#'     in the data.
#' @param splinesobj list of spline components in the MTRs for treated
#'     and control groups. Spline terms are extracted using
#'     \code{\link{removeSplines}}. This object is supposed to be a
#'     dictionary of splines, containing the original calls of each
#'     spline in the MTRs, their specifications, and the index used
#'     for renaming each component.
#' @param environments a list containing the environments of the MTR
#'     formulas, the IV-like formulas, and the propensity score
#'     formulas. If a formulas is not provided, and thus no
#'     environment can be found, then the parent.frame() is assigned
#'     by default.
#'
#' @inheritParams ivmte
#'
#' @return Returns a list of results from throughout the estimation
#'     procedure. This includes all IV-like estimands; the propensity
#'     score model; bounds on the treatment effect; the estimated
#'     expectations of each term in the MTRs; the components and
#'     results of the LP problem.
# @export
ivmteEstimate <- function(data, target, late.Z, late.from, late.to,
                          late.X, eval.X, genlate.lb, genlate.ub,
                          target.weight0, target.weight1,
                          target.knots0 = NULL, target.knots1 = NULL,
                          m0, m1, uname = u, m1.ub, m0.ub, m1.lb,
                          m0.lb, mte.ub, mte.lb, m0.dec, m0.inc,
                          m1.dec, m1.inc, mte.dec, mte.inc, ivlike,
                          components, subset, propensity,
                          link = "logit", treat, lpsolver,
                          lpsolver.options, lpsolver.presolve,
                          lpsolver.options.criterion, lpsolver.options.bounds,
                          criterion.tol = 0, initgrid.nx = 20,
                          initgrid.nu = 20, audit.nx = 2500,
                          audit.nu = 25, audit.add = 100,
                          audit.max = 25, audit.tol, audit.grid = NULL,
                          point = FALSE,
                          point.eyeweight = FALSE,
                          point.center = NULL, point.redundant = NULL,
                          count.moments = TRUE,
                          orig.sset = NULL,
                          orig.criterion = NULL, vars_y,
                          vars_mtr, terms_mtr0, terms_mtr1, vars_data,
                          splinesobj, noisy = TRUE,
                          smallreturnlist = FALSE, seed = 12345,
                          debug = FALSE, environments) {
    call <- match.call(expand.dots = FALSE)
    if (classFormula(ivlike)) ivlike <- c(ivlike)
    ## Character arguments will be converted to lowercase
    if (hasArg(lpsolver)) lpsolver <- tolower(lpsolver)
    if (hasArg(target))   target   <- tolower(target)
    if (hasArg(link))     link     <- tolower(link)
    if (hasArg(ci.type))  ci.type  <- tolower(ci.type)
    ## Convert uname into a string
    uname <- deparse(substitute(uname))
    uname <- gsub("~", "", uname)
    uname <- gsub("\\\"", "", uname)
    if (noisy == TRUE && hasArg(lpsolver)) {
        if (lpsolver == "gurobi") cat("\nLP solver: Gurobi ('gurobi')\n\n")
        if (lpsolver == "cplexapi") cat("\nLP solver: CPLEX ('cplexAPI')\n\n")
        if (lpsolver == "lpsolveapi") {
            cat("\nLP solver: lp_solve ('lpSolveAPI')\n\n")
            warning(gsub("\\s+", " ",
                     "The R package 'lpSolveAPI' interfaces with 'lp_solve',
                      which is outdated and potentially unreliable. It is
                      recommended to use commercial solvers
                      Gurobi (lpsolver = 'gurobi')
                      or CPLEX (lpsolver = 'cplexAPI') instead.
                      Free academic licenses can be obtained for these
                      commercial solvers."),
                "\n", call. = FALSE, immediate. = TRUE)
        }
    }

    ##---------------------------
    ## 1. Obtain propensity scores
    ##---------------------------
    if (noisy == TRUE) {
        cat("Obtaining propensity scores...\n")
    }
    ## Estimate propensity scores
    pcall <- modcall(call,
                     newcall = propensity,
                     keepargs = c("link", "late.Z", "late.X"),
                     dropargs = "propensity",
                     newargs = list(data = quote(data),
                                    formula = propensity,
                                    env = environments$propensity))
    pmodel <- eval(pcall)

    ##---------------------------
    ## 2. Generate target moments/gamma terms
    ##---------------------------
    if (noisy == TRUE) {
        cat("\nGenerating target moments...\n")
    }
    ## Parse polynomials
    if (!is.null(m0)) {
        m0call <- modcall(call,
                          newcall = polyparse,
                          keepargs = c("uname"),
                          newargs = list(formula = m0,
                                         data = quote(data),
                                         env = quote(environments$m0)))
        pm0 <- eval(as.call(m0call))
    } else {
        pm0 <- NULL
    }
    if (!is.null(m1)) {
        m1call <- modcall(call,
                          newcall = polyparse,
                          keepargs = c("uname"),
                          newargs = list(formula = m1,
                                         data = quote(data),
                                         env = quote(environments$m1)))
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
                                              "genlate.ub", "splinesobj",
                                              "point", "noisy"),
                                 dropargs = "data",
                                 newargs = list(data = quote(data),
                                                pmodobj = quote(pmodel),
                                                pm0 = quote(pm0),
                                                pm1 = quote(pm1)))
    } else {
        gentargetcall <- modcall(call,
                                 newcall = genTarget,
                                 keepargs = c("treat", "m1", "m0",
                                              "target.weight0",
                                              "target.weight1",
                                              "target.knots0", "target.knots1",
                                              "splinesobj",
                                              "point", "noisy"),
                                 dropargs = "data",
                                 newargs = list(data = quote(data),
                                                pmodobj = quote(pmodel),
                                                pm0 = quote(pm0),
                                                pm1 = quote(pm1)))
    }
    targetGammas <- eval(gentargetcall)
    rm(gentargetcall)
    gstar0 <- targetGammas$gstar0
    gstar1 <- targetGammas$gstar1

    ##---------------------------
    ## 3. Generate moments/gamma terms for IV-like estimands
    ##---------------------------
    if (noisy == TRUE) {
        cat("\nGenerating IV-like moments...\n")
    }
    sset  <- list() ## Contains all IV-like estimates and their
                    ## corresponding moments/gammas
    scount <- 1     ## counter for S-set constraints
    subsetIndexList <- list()
    if (classList(ivlike)) {
        ## Construct `sset' object when multiple IV-like
        ## specifications are provided
        ivlikeCounter <- 1
        ivlikeD <- NULL
        for (i in 1:length(ivlike)) {
            sformula   <- ivlike[[i]]
            if (all(length(Formula::as.Formula(sformula)) == c(1, 1))) {
                if (treat %in% all.vars(Formula::as.Formula(sformula)[[3]])) {
                    ivlikeD <- c(ivlikeD, TRUE)
                } else {
                    ivlikeD <- c(ivlikeD, FALSE)
                }
            }
            if (all(length(Formula::as.Formula(sformula)) == c(1, 2))) {
                if (treat %in%
                    all.vars(Formula::as.Formula(sformula)[[3]][[2]])) {
                    ivlikeD <- c(ivlikeD, TRUE)
                } else {
                    ivlikeD <- c(ivlikeD, FALSE)
                }
            }
            environment(sformula) <- environments$ivlike
            scomponent <- components[[i]]
            if (subset[[i]] == "") {
                ssubset <- replicate(nrow(data), TRUE)
            } else {
                ssubset <- subset[[i]]
            }
            ## Obtain coefficient estimates and S-weights
            ## corresponding to the IV-like estimands
            sest  <- ivEstimate(formula = sformula,
                                data = data[eval(substitute(ssubset), data), ],
                                components = scomponent,
                                treat = treat,
                                list = TRUE,
                                order = ivlikeCounter)
            ## Generate moments (gammas) corresponding to IV-like
            ## estimands
            subset_index <- rownames(data[eval(substitute(ssubset), data), ])
            if (length(subset_index) == nrow(data)) {
                subsetIndexList[[ivlikeCounter]] <- NA
            } else {
                subsetIndexList[[ivlikeCounter]] <- as.integer(subset_index)
            }
            ncomponents <- sum(!is.na(sest$betas))
            if (point) {
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
                                  noisy = noisy,
                                  ivn = ivlikeCounter,
                                  redundant = point.redundant)
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
                                  means = TRUE,
                                  yvar = vars_y,
                                  dvar = treat,
                                  noisy = noisy,
                                  ivn = ivlikeCounter,
                                  redundant = point.redundant)
            }
            ## Update set of moments (gammas)
            sset <- setobj$sset
            scount <- setobj$scount
            rm(setobj)
            ivlikeCounter <- ivlikeCounter + 1
        }
    } else {
        stop(gsub("\\s+", " ",
                  "'ivlike' argument must either be a formula or a vector of
                  formulas."),
             call. = FALSE)
    }
    rm(sest, subset_index)
    if (exists("pm0")) rm(pm0)
    if (exists("pm1")) rm(pm1)
    if (smallreturnlist) pmodel <- pmodel$model
    if (count.moments) {
        wmat <- NULL
        for (s in 1:length(sset)) {
            if (!is.null(subsetIndexList)) {
                wmatTmp <- rep(0,  times = 2 * nrow(data))
                if (!is.integer(subsetIndexList[[sset[[s]]$ivspec]])) {
                    wmatTmp <- c(sset[[s]]$w0$multiplier,
                                 sset[[s]]$w1$multiplier)
                } else {
                    wmatTmp[c(subsetIndexList[[sset[[s]]$ivspec]],
                              nrow(data) +
                              subsetIndexList[[sset[[s]]$ivspec]])] <-
                        c(sset[[s]]$w0$multiplier, sset[[s]]$w1$multiplier)
                }
            } else {
                wmatTmp <- c(sset[[s]]$w0$multiplier,
                             sset[[s]]$w1$multiplier)
            }
            wmat <- cbind(wmat, wmatTmp)
            rm(wmatTmp)
        }
        ## Check number of linearly independent moments
        nIndepMoments <- qr(wmat)$rank
        if (noisy == TRUE) cat("    Independent moments:", nIndepMoments, "\n")
        if (nIndepMoments < length(sset) && !all(ivlikeD)) {
            warning(gsub("\\s+", " ",
                         paste0("The following IV-like specifications do not
                            include the treatment variable: ",
                            paste(which(!ivlikeD), collapse = ", "),
                            ". This may result in fewer
                            independent moment conditions than expected.")),
                    call. = FALSE)
        }
        rm(wmat)
    } else {
        nIndepMoments <- NULL
    }
    if (!is.null(point.redundant)) point.redundant <- 0
    ## If bootstrapping, check that length of sset is equivalent in
    ## length to that of the original sset if bootstrapping
    if (!is.null(orig.sset)) {
        if (length(sset) != length(orig.sset)) {
            return("collinearity causing S-set to differ from original sample")
        }
    }
    ## Prepare GMM estimate estimate if `point' agument is set to TRUE
    splinesCheck <- !(all(is.null(splinesobj[[1]]$splineslist)) &&
        all(is.null(splinesobj[[2]]$splineslist)))
    if (point == TRUE) {
        ## Obtain GMM estimate
        gmmResult <- gmmEstimate(sset = sset,
                                 gstar0 = gstar0,
                                 gstar1 = gstar1,
                                 center = point.center,
                                 identity = point.eyeweight,
                                 redundant = point.redundant,
                                 subsetList = subsetIndexList,
                                 n = nrow(data),
                                 nMoments = nIndepMoments,
                                 splines = splinesCheck,
                                 noisy = noisy)
        if (!smallreturnlist) {
            return(list(sset  = sset,
                        gstar = list(g0 = colMeans(gstar0),
                                     g1 = colMeans(gstar1)),
                        propensity = pmodel,
                        pointestimate = gmmResult$pointestimate,
                        moments = gmmResult$moments,
                        redundant = gmmResult$redundant,
                        jtest = gmmResult$jtest,
                        mtr.coef = gmmResult$coef))
        } else {
            sset <- lapply(sset, function(x) {
                output <- list()
                output$ivspec <- x$ivspec
                output$beta <- x$beta
                output$g0 <- colMeans(x$g0)
                output$g1 <- colMeans(x$g1)
                return(output)
            })
            output <- list(sset = sset,
                           gstar = list(g0 = colMeans(gstar0),
                                        g1 = colMeans(gstar1)),
                           pointestimate = gmmResult$pointestimate,
                           moments = gmmResult$moments,
                           redundant = gmmResult$redundant,
                           jtest = gmmResult$jtest,
                           mtr.coef = gmmResult$coef)
            if (all(class(pmodel) != "NULL")) {
                output$propensity.coef <- pmodel$coef
            }
            return(output)
        }
    }

    ##---------------------------
    ## 4. Define constraint matrices using the audit
    ##---------------------------
    if (noisy == TRUE) {
        cat("\nPerforming audit procedure...\n")
    }
    audit.args <- c("uname", "vars_data",
                    "initgrid.nu", "initgrid.nx",
                    "audit.nx", "audit.nu", "audit.add",
                    "audit.max", "audit.tol",
                    "audit.grid",
                    "m1.ub", "m0.ub",
                    "m1.lb", "m0.lb",
                    "mte.ub", "mte.lb", "m0.dec",
                    "m0.inc", "m1.dec", "m1.inc", "mte.dec",
                    "mte.inc", "lpsolver.options", "lpsolver.presolve",
                    "lpsolver.options.criterion", "lpsolver.options.bounds",
                    "criterion.tol",
                    "orig.sset", "orig.criterion",
                    "smallreturnlist", "noisy", "seed", "debug")
    audit_call <- modcall(call,
                          newcall = audit,
                          keepargs = audit.args,
                          newargs = list(data = quote(data),
                                         m0 = quote(m0),
                                         m1 = quote(m1),
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
    autoExpand <- 0
    autoExpandMax <- 3
    newGrid.nu <- initgrid.nu
    newGrid.nx <- initgrid.nx
    while(autoExpand <= autoExpandMax) {
        audit <- eval(audit_call)
        if (is.null(audit$error)) {
            autoExpand <- Inf
        }
        if (!is.null(audit$error) &&
            audit$error == "Failure to maximize/minimize.") {
            if (initgrid.nx == audit.nx && initgrid.nu == audit.nu) {
                stop(paste0(gsub("\\s+", " ",
                                 "The LP problem is unbounded,
                                  and automatic grid expansion is not
                                  possible.
                                  Consider imposing additional shape
                                  constraints, or increasing the size of the
                                  initial grid  for the audit.
                                  In order to increase the size of the
                                  initial grid, it will be
                                  necessary to increase the size of the audit
                                  grid.
                                  This is because the initial grid must be a
                                  subset of the audit grid, and this constraint
                                  is currently binding.
                                  In order to allow for automatic grid
                                  expansion, make sure audit.nx > initgrid.nx
                                  or audit.nu > initgrid.nu.")),
                     call. = FALSE)
            }
            autoExpand <- autoExpand + 1
            newGrid.nu <- min(ceiling(newGrid.nu * 1.5), audit.nu)
            newGrid.nx <- min(ceiling(newGrid.nx * 1.5), audit.nx)
            cat("\n    Restarting audit with expanded initial grid.\n")
            cat(paste0("    New settings: initgrid.nx = ",  newGrid.nx,
                       ", initgrid.nu = ", newGrid.nu, "\n"))
            audit_call <-
                modcall(audit_call,
                        dropargs = c("initgrid.nu", "initgrid.nx",
                                     "audit.grid"),
                        newargs = list(initgrid.nu = newGrid.nu,
                                       initgrid.nx = newGrid.nx,
                                       audit.grid = audit$audit.grid))
            if (newGrid.nu == audit.nu && newGrid.nx == audit.nx) {
                autoExpand <- autoExpandMax
            }
            if (autoExpand == autoExpandMax) {
                autoExpand <- Inf
                audit <- eval(audit_call)
            }
        }
    }
    rm(audit_call)
    if (!is.null(audit$error) && autoExpand > autoExpandMax) {
        stop(paste0(gsub("\\s+", " ",
                         "Automatic grid expansion limit reached.
                         The LP problem is still unbounded.
                         Either impose additional shape constraints,
                         or increase the size of the initial grid
                         for the audit. Since the initial grid must
                         be a subset of the audit grid, it may be
                         necessary to increase the size of the audit
                         grid also.
                         The most recent options after automatic expansion
                         were:"),
                    "\ninitgrid.nx = ", newGrid.nx,
                    "\ninitgrid.nu = ", newGrid.nu,
                    "\naudit.nx = ", audit.nx,
                    "\naudit.nu = ", audit.nu),
             call. = FALSE)
    }
    if (noisy) {
        cat("Bounds on the target parameter: [",
            fmtResult(audit$min), ", ", fmtResult(audit$max), "]\n\n", sep = "")
        ## if (any(audit$lpresult$modelstats[, 3] > 6)) {
        ##     bMessage <- "The following sets of coefficients defining the
        ##             LP problem exhibit ranges exceeding 6 orders of magnitude: "
        ##     if (audit$lpresult$modelstats[1, 3] > 6) {
        ##         bMessage <- paste(bMessage, "constraint matrix")
        ##     }
        ##     if (audit$lpresult$modelstats[2, 3] > 6) {
        ##         bMessage <- paste(bMessage, "RHS vector (IV-like coefficients)")
        ##     }
        ##     if (audit$lpresult$modelstats[3, 3] > 6) {
        ##         bMessage <- paste(bMessage, "objective vector (gamma moments)")
        ##     }
        ##     bMessage <- paste0(bMessage, ". Large ranges in the coefficients
        ##                                  increase computational burden, and can
        ##                                  potentially lead to infeasibility.")
        ##     warning(gsub("\\s+", " ", bMessage),
        ##             call. = FALSE, immediate. = TRUE)
        ## }
    }
    ## include additional output material
    if (lpsolver == "gurobi") lpsolver <- "Gurobi ('gurobi')"
    if (lpsolver == "lpsolveapi") lpsolver <- "lp_solve ('lpSolveAPI')"
    if (lpsolver == "cplexapi") lpsolver <- "CPLEX ('cplexAPI')"
    if (!smallreturnlist) {
        output <- list(sset  = sset,
                       gstar = list(g0 = gstar0,
                                    g1 = gstar1,
                                 n = targetGammas$n),
                       gstar.weights = list(w0 = targetGammas$w0,
                                         w1 = targetGammas$w1),
                       gstar.coef = list(min.g0 = audit$lpresult$ming0,
                                         max.g0 = audit$lpresult$maxg0,
                                         min.g1 = audit$lpresult$ming1,
                                         max.g1 = audit$lpresult$maxg1),
                       propensity = pmodel,
                       bounds = c(audit$min, audit$max),
                       lpresult =  audit$lpresult,
                       lpsolver = lpsolver,
                       indep.moments = nIndepMoments,
                       audit.grid = list(audit.x =
                                             audit$gridobj$audit.grid$support,
                                         audit.u =
                                             audit$gridobj$audit.grid$uvec,
                                         audit.noX =
                                             audit$gridobj$audit.grid$noX,
                                         violations = audit$gridobj$violations),
                       audit.count = audit$auditcount,
                       audit.criterion = audit$minobseq,
                       splinesdict = list(m0 = splinesobj[[1]]$splinesdict,
                                          m1 = splinesobj[[2]]$splinesdict))
    } else {
        sset <- lapply(sset, function(x) {
            x[c("ivspec", "beta", "g0", "g1")]
        })
        output <- list(sset  = sset,
                       gstar = list(g0 = gstar0,
                                    g1 = gstar1),
                       gstarcoef = list(ming0 = audit$ming0,
                                        maxg0 = audit$maxg0,
                                        ming1 = audit$ming1,
                                        maxg1 = audit$maxg1),
                       bounds = c(audit$min, audit$max),
                       lpsolver = lpsolver,
                       indep.moments = nIndepMoments,
                       audit.count = audit$auditcount,
                       audit.criterion = audit$minobseq,
                       splinesdict = list(m0 = splinesobj[[1]]$splinesdict,
                                          m1 = splinesobj[[2]]$splinesdict))
        if (all(class(pmodel$model) != "NULL")) {
            output$propensity.coef <- pmodel$model$coef
        }
    }
    if (!is.null(audit$spectest)) output$specification.test <- audit$spectest
    return(output)
}


#' Generating LP moments for IV-like estimands
#'
#' This function takes in the IV estimate and its IV-like
#' specification, and generates a list containing the corresponding
#' point estimate, and the corresponding moments (gammas) that will
#' enter into the constraint matrix of the LP problem.
#'
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
#'
#' @inheritParams ivmteEstimate
#'
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
#' dtm <- ivmte:::gendistMosquito()
#'
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
genTarget <- function(treat, m0, m1, target,
                      target.weight0, target.weight1,
                      target.knots0, target.knots1,
                      late.Z, late.from, late.to, late.X,
                      eval.X, genlate.lb, genlate.ub,
                      data, splinesobj, pmodobj, pm0, pm1,
                      point = FALSE, noisy = TRUE) {
    if (!hasArg(late.X)) {
        late.X <- NULL
        eval.X <- NULL
        lateRows <- NULL
    }
    if (hasArg(target)) target   <- tolower(target)
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
            if (!is.null(late.X)) {
                ## Create index to restrict data to those with the
                ## correct value of X
                condX <- mapply(function(a, b) paste(a, "==", b),
                                late.X, eval.X)
                condX <- paste(condX, collapse = " & ")
                lateRows <- eval(parse(text = condX), data)
                if (sum(lateRows) == 0) {
                    stop(gsub("\\s+", " ",
                              "no observations with the values specified in
                               'eval.X'."), call. = FALSE)
                }
            }
            if (!is.null(m0)) pm0$polymat <- pm0$polymat
            if (!is.null(m1)) pm1$polymat <- pm1$polymat
            w1 <- wlate1(data, late.from, late.to, late.Z,
                         pmodobj$model, late.X, eval.X)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else if (target == "genlate") {
            if (!is.null(late.X)) {
                ## Create index to restrict data to those with the
                ## correct value of X
                condX <- mapply(function(a, b) paste(a, "==", b),
                                late.X, eval.X)
                condX <- paste(condX, collapse = " & ")
                lateRows <- eval(parse(text = condX), data)
                if (sum(lateRows) == 0) {
                    stop(gsub("\\s+", " ",
                              "no observations with the values specified in
                               'eval.X'."), call. = FALSE)
                }
            }
            if (!is.null(m0)) pm0$polymat <- pm0$polymat
            if (!is.null(m1)) pm1$polymat <- pm1$polymat
            w1 <- wgenlate1(data, genlate.lb, genlate.ub)
            w0 <- w1
            w0$mp <- -1 * w0$mp
        } else {
            stop("Unrecognized target parameter.",
                 call. = FALSE)
        }
        ## Integrate m0 and m1 functions
        if (!is.null(m0)) {
            if (noisy == TRUE) {
                cat("    Integrating terms for control group...\n")
            }
            if (point == FALSE) {
                gstar0 <- genGamma(monomials = pm0,
                                   lb = w0$lb,
                                   ub = w0$ub,
                                   multiplier = w0$mp,
                                   late.rows = lateRows)
            } else {
                gstar0 <- genGamma(monomials = pm0,
                                   lb = w0$lb,
                                   ub = w0$ub,
                                   multiplier = w0$mp,
                                   means = FALSE,
                                   late.rows = lateRows)
            }
        } else {
            gstar0 <- NULL
        }
        if (!is.null(m1)) {
            if (noisy == TRUE) {
                cat("    Integrating terms for treated group...\n")
            }
            if (point == FALSE) {
                gstar1 <- genGamma(monomials = pm1,
                                   lb = w1$lb,
                                   ub = w1$ub,
                                   multiplier = w1$mp,
                                   late.rows = lateRows)
            } else {
                gstar1 <- genGamma(monomials = pm1,
                                   lb = w1$lb,
                                   ub = w1$ub,
                                   multiplier = w1$mp,
                                   means = FALSE,
                                   late.rows = lateRows)
            }
        } else {
            gstar1 <- NULL
        }
        if (point == FALSE) {
            gstarSplineObj0 <- genGammaSplines(splinesobj = splinesobj[[1]],
                                               data = data,
                                               lb = w0$lb,
                                               ub = w0$ub,
                                               multiplier = w0$mp,
                                               d = 0,
                                               late.rows = lateRows)
            gstarSpline0 <- gstarSplineObj0$gamma
            gstarSplineObj1 <- genGammaSplines(splinesobj = splinesobj[[2]],
                                               data = data,
                                               lb = w1$lb,
                                               ub = w1$ub,
                                               multiplier = w1$mp,
                                               d = 1,
                                               late.rows = lateRows)
            gstarSpline1 <- gstarSplineObj1$gamma
        } else {
            gstarSplineObj0 <- genGammaSplines(splinesobj = splinesobj[[1]],
                                               data = data,
                                               lb = w0$lb,
                                               ub = w0$ub,
                                               multiplier = w0$mp,
                                               d = 0,
                                               means = FALSE,
                                               late.rows = lateRows)
            gstarSpline0 <- gstarSplineObj0$gamma
            gstarSplineObj1 <- genGammaSplines(splinesobj = splinesobj[[2]],
                                               data = data,
                                               lb = w1$lb,
                                               ub = w1$ub,
                                               multiplier = w1$mp,
                                               d = 1,
                                               means = FALSE,
                                               late.rows = lateRows)
            gstarSpline1 <- gstarSplineObj1$gamma
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
                        cat("    ",
                            gsub("\\s+", " ",
                                 "Integrating non-spline
                                  terms for control group..."),
                            "\n", sep = "")
                    } else {
                        cat("    ",
                            gsub("\\s+", " ",
                                 "Integrating non-spline
                                  terms for treated group..."),
                            "\n", sep = "")
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
                    cat("    Integrating spline terms for control group...\n")
                    } else {
                    cat("    Integrating spline terms for treated group...\n")
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
                            genGammaSplines(splinesobj = noSplineMtr,
                                            data = data,
                                            lb = lb,
                                            ub = ub,
                                            multiplier = weights,
                                            d = d)$gamma
                    } else {
                        gammaSplines <- gammaSplines +
                            genGammaSplines(splinesobj = noSplineMtr,
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
                   gstar1 = gstar1)
    if (hasArg(target)) {
        output$w1 <- w1
        output$w0 <- w0
    }
    output$n <- nrow(data)
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
#' @param redundant vector of integers indicating which components in
#'     the S-set are redundant.
#' @return A list containing the point estimate for the IV regression,
#'     and the expectation of each monomial term in the MTR.
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#'
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
                    yvar, dvar, noisy = TRUE, ivn = NULL,
                    redundant = NULL) {
    if (!hasArg(subset_index)) {
        subset_index <- NULL
        n <- nrow(data)
    } else {
        n <- length(subset_index)
    }
    for (j in 1:ncomponents) {
        if (! scount %in% redundant) {
            if (noisy == TRUE) {
                cat("    Moment ", scount, "...\n", sep = "")
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
                gsSpline0 <- genGammaSplines(splinesobj = splinesobj[[1]],
                                             data = data,
                                             lb = pmodobj,
                                             ub = 1,
                                             multiplier = sest$sw0[, j],
                                             subset = subset_index,
                                             d = 0)$gamma
                gsSpline1 <- genGammaSplines(splinesobj = splinesobj[[2]],
                                             data = data,
                                             lb = 0,
                                             ub = pmodobj,
                                             multiplier = sest$sw1[, j],
                                             subset = subset_index,
                                             d = 1)$gamma
            } else {
                gsSpline0 <- genGammaSplines(splinesobj = splinesobj[[1]],
                                             data = data,
                                             lb = pmodobj,
                                             ub = 1,
                                             multiplier = sest$sw0[, j],
                                             subset = subset_index,
                                             d = 0,
                                             means = FALSE)$gamma
                gsSpline1 <- genGammaSplines(splinesobj = splinesobj[[2]],
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
                                                    w1 = sweight1,
                                                    n = n)
            } else {
                ## Now generate the vectors for Y * S(D, Z).
                if (!is.null(subset_index)) {
                    newsubset <- subset_index
                } else {
                    newsubset <- seq(1, nrow(data))
                }
                yvec <- as.vector(data[newsubset, yvar])
                dvec <- as.vector(data[newsubset, dvar])
                yvec <- yvec * (sest$sw1[, j] * dvec + sest$sw0[, j] *
                                (1 - dvec))
                names(yvec) <- newsubset
                sset[[paste0("s", scount)]] <- list(ivspec = ivn,
                                                    beta = sest$beta[j],
                                                    g0 = cbind(gs0, gsSpline0),
                                                    g1 = cbind(gs1, gsSpline1),
                                                    ys = yvec,
                                                    w0 = sweight0,
                                                    w1 = sweight1,
                                                    n = n)
            }
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
#' @param center numeric, the GMM moment equations from the original
#'     sample. When bootstrapping, the solution to the point
#'     identified case obtained from the original sample can be passed
#'     through this argument to recenter the bootstrap distribution of
#'     the J-statistic.
#' @param subsetList list of subset indexes, one for each IV-like
#'     specification.
#' @param n number of observations in the data. This option is only
#'     used when subsetting is involved.
#' @param redundant vector of integers indicating which components in
#'     the S-set are redundant.
#' @param identity boolean, default set to \code{FALSE}. Set to
#'     \code{TRUE} if GMM point estimate should use the identity
#'     weighting matrix (i.e. one-step GMM).
#' @param nMoments number of linearly independent moments. This option
#'     is used to determine the cause of underidentified cases.
#' @param splines boolean, set to \code{TRUE} if the MTRs involve
#'     splines. This option is used to determine the cause of
#'     underidentified cases.
#' @param noisy boolean, default set to \code{TRUE}. If \code{TRUE},
#'     then messages are provided throughout the estimation
#'     procedure. Set to \code{FALSE} to suppress all messages,
#'     e.g. when performing the bootstrap.
#' @return a list containing the point estimate of the treatment
#'     effects, and the MTR coefficient estimates. The moment
#'     conditions evaluated at the solution are also returned, along
#'     with the J-test results. However, if the option \code{center}
#'     is passed, then the moment conditions and J-test are centered.
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#'
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
gmmEstimate <- function(sset, gstar0, gstar1, center = NULL,
                        subsetList = NULL, n = NULL, redundant = NULL,
                        identity = FALSE, nMoments, splines,
                        noisy = TRUE) {
    gn0 <- ncol(gstar0)
    gn1 <- ncol(gstar1)
    if ((gn0 + gn1) > length(sset)) {
        stop(gsub("\\s+", " ",
                  paste0("GMM system is underidentified: there are ",
                         (gn0 + gn1),
                         " unknown MTR coefficients and ",
                         length(sset),
                         " moment conditions defined by IV-like
                         specifications. Either expand the number of
                         moments, or adjust m0 and m1.")),
             call. = FALSE)
    }
    ## Perform first step of GMM
    mm <- momentMatrix(sset, gn0, gn1, subsetList, n)
    mlist <- (seq(length(sset)) - 1) * (gn0 + gn1 + 1) + 1
    altmm <- matrix(t(mm[, -mlist]), byrow = TRUE, ncol = (gn0 + gn1))
    altmean <- Matrix::Matrix(t(rep(1, nrow(mm)) %x% diag(length(sset))),
                              sparse = TRUE)
    xmat <- matrix(altmean %*% altmm / nrow(mm), ncol = (gn0 + gn1))
    if (qr(xmat)$rank < (gn0 + gn1)) {
        if (hasArg(nMoments) && hasArg(splines)) {
            if (nMoments > qr(xmat)$rank && splines) {
                stop(gsub("\\s+", " ",
                          paste0("GMM system is underidentified: there are ",
                          (gn0 + gn1),
                          " unknowns, but the Gamma matrix has rank ",
                          qr(xmat)$rank,
                          ". The previous count of ",
                          nMoments,
                          " independent moments becomes
                          unreliable when m0 or m1 includes splines.
                          Either adjust the IV-like specifications,
                          or adjust m0 and m1.")),
                     call. = FALSE)
            }
        }
        stop(gsub("\\s+", " ",
                  paste0("GMM system is underidentified: there are ",
                  (gn0 + gn1),
                  " unknown MTR coefficients and ",
                  qr(xmat)$rank,
                  " linearly independent moment conditions defined by the
                    IV-like specifications. Either adjust the
                    IV-like specifications, or adjust m0 and m1.")),
             call. = FALSE)
    }
    ymat <- colMeans(mm[, mlist])
    ## Address collinear moments
    if (is.null(redundant)) {
        colDrop <- NULL
        theta <- rnorm(gn0 + gn1) ## Random theta, to test for collinearity
        tmpOmegaMat <- c(t(mm[, mlist]))
        tmpOmegaMat <- tmpOmegaMat - altmm %*% theta
        tmpOmegaMat <- matrix(tmpOmegaMat, byrow = TRUE, ncol = length(sset))
        tmpOmegaMat <- t(tmpOmegaMat) %*% tmpOmegaMat / nrow(mm)
        rankCheck <- eigen(tmpOmegaMat)
        if (any(abs(rankCheck$values) < 1e-08)) {
            colDict <- list()
            colDrop <- seq(ncol(tmpOmegaMat))
            colnames(tmpOmegaMat) <- seq(ncol(tmpOmegaMat))
            while (any(abs(rankCheck$values) < 1e-08)) {
                colPos <- which(abs(rankCheck$values) < 1e-08)[1]
                colVec <- rankCheck$vectors[, colPos]
                colSeq <- which(abs(colVec) > 1e-08)
                colDropIndex <- max(colSeq)
                tmpOmegaMat <- tmpOmegaMat[-colDropIndex, -colDropIndex]
                rankCheck <- eigen(tmpOmegaMat)
            }
            colDrop <- which(! colDrop %in% colnames(tmpOmegaMat))
        }
    } else {
        if (length(redundant) == 1 && redundant == 0) {
            colDrop <- NULL
        } else {
            colDrop <- redundant
        }
    }
    if (!is.null(colDrop)) {
        colDrop <- sort(colDrop)
        xmat <- xmat[-colDrop, ]
        ymat <- ymat[-colDrop]
        colDropStr <- sapply(colDrop, FUN = function(x) {
            c(sset[[x]]$ivspec, names(sset[[x]]$beta))
        })
        if (is.null(redundant)) {
            dropMessage <- paste("IV-like Spec.     Component\n",
                                 "------------------------------------------\n",
                                 collapse = "")
            for (i in 1:ncol(colDropStr)) {
                dropMessage <-
                    paste(dropMessage, colDropStr[1, i],
                          paste(rep(" ", 16 - length(colDropStr[1, i])),
                                collapse = ""),
                          colDropStr[2, i], "\n",
                          collapse = "")
            }
            warning(paste(gsub("\\s+", " ",
                               "The following components have been dropped
                                due to collinearity:"),
                          "\n",
                          dropMessage,
                          collapse = ""),
                    call. = FALSE,
                    immediate. = FALSE)
        }
    }
    else {
        if (length(redundant) == 1 && redundant == 0) {
            colDrop <- NULL
        } else {
            colDrop <- redundant
        }
    }
    ## Perform first stage GMM
    if (is.null(center)) {
        theta <- solve(t(xmat) %*% xmat) %*% t(xmat) %*% ymat
    } else {
        theta <- solve(t(xmat) %*% xmat) %*%
            (t(xmat) %*% ymat - t(xmat) %*% center)
    }
    ## Perform second step of GMM
    if (!identity) {
        alty <- c(t(mm[, mlist]))
        omegaMat <- alty - altmm %*% theta
        omegaMat <- matrix(omegaMat, byrow = TRUE, ncol = length(sset))
        omegaMat <- t(omegaMat) %*% omegaMat / nrow(mm)
        if (!is.null(colDrop)) omegaMat <- omegaMat[-colDrop, -colDrop]
        omegaMat <- solve(omegaMat)
        if (is.null(center)) {
            theta <- solve(t(xmat) %*% omegaMat %*% xmat) %*%
                t(xmat) %*% omegaMat %*% ymat
        } else {
            theta <- solve(t(xmat) %*% omegaMat %*% xmat) %*%
                (t(xmat) %*% omegaMat %*% ymat -
                 t(xmat) %*% omegaMat %*% center)
        }
    }
    moments <- ymat - xmat %*% theta
    ## Perform J test
    if ((length(sset) - length(colDrop)) > gn0 + gn1) {
        if (identity) {
            if (is.null(center)) {
                jtest <- nrow(mm) * t(moments) %*% moments
            } else {
                jtest <- nrow(mm) * t(moments - center) %*%
                    (moments - center)
            }
        } else {
            if (is.null(center)) {
                jtest <- nrow(mm) * t(moments) %*% omegaMat %*% moments
            } else {
                jtest <- nrow(mm) * t(moments - center) %*%
                    omegaMat %*% (moments - center)
            }
        }
        if (is.null(redundant)) {
            chiDf <- length(sset) - length(colDrop) - gn0 - gn1
        } else {
            chiDf <- length(sset) - gn0 - gn1
        }
        jtest <- c(jtest,
                   1 - pchisq(jtest,
                              df = chiDf),
                   chiDf)
        names(jtest) <- c("J-statistic",
                          "p-value (ignoring first step)", "df")
    } else {
        jtest <- NULL
    }
    names(theta) <- c(paste0("m0.", colnames(gstar0)),
                      paste0("m1.", colnames(gstar1)))
    pointestimate <- sum(c(colMeans(gstar0), colMeans(gstar1)) * theta)
    if (noisy == TRUE) {
        cat("\nPoint estimate of the target parameter: ",
            pointestimate, "\n\n", sep = "")
    }
    if (is.null(colDrop)) colDrop <- 0
    return(list(pointestimate = as.numeric(pointestimate),
                coef = theta,
                moments = moments,
                jtest = jtest,
                redundant = colDrop))
}

#' Construct pre-meaned moment matrix
#'
#' This function constructs the matrix to be fed into the GMM
#' estimator to construct the moment conditions.
#'
#' @param sset a list of lists constructed from the function
#'     \link{genSSet}. Each inner list should include a coefficient
#'     corresponding to a term in an IV specification, a matrix of the
#'     estimates of the gamma moments conditional on (X, Z) for d = 0,
#'     and a matrix of the estimates of the gamma moments conditional
#'     on (X, Z) for d = 1. The column means of the last two matrices
#'     is what is used to generate the gamma moments.
#' @param gn0 integer, number of terms in the MTR for control group.
#' @param gn1 integer, number of terms in the MTR for treated group.
#' @param subsetList list of subset indexes, one for each IV-like
#'     specification.
#' @param n number of observations in the data. This option is only
#'     used when subsets are involved.
#' @return matrix whose column means can be used to carry out the GMM
#'     estimation.
momentMatrix <- function(sset, gn0, gn1, subsetList = NULL, n = NULL) {
    momentMatrix <- NULL
    momentNames <- NULL
    for (s in 1:length(sset)) {
        if (!is.null(subsetList)) {
            momentMatrixTmp <- matrix(0, nrow = n, ncol = 1 + gn0 + gn1)
            if (!is.integer(subsetList[[sset[[s]]$ivspec]])) {
                momentMatrixTmp <- cbind(sset[[s]]$ys,
                                         sset[[s]]$g0,
                                         sset[[s]]$g1)
            } else {
                momentMatrixTmp[subsetList[[sset[[s]]$ivspec]], ] <-
                    cbind(sset[[s]]$ys , sset[[s]]$g0, sset[[s]]$g1)
            }
        } else {
            momentMatrixTmp <- cbind(sset[[s]]$ys, sset[[s]]$g0, sset[[s]]$g1)
        }
        momentMatrix <- cbind(momentMatrix,
                              momentMatrixTmp)
        momentNames <- c(momentNames,
                         paste0("s", s, "y"),
                         paste0("s", s, "g0", seq(1, gn0)),
                         paste0("s", s, "g1", seq(1, gn1)))
    }
    colnames(momentMatrix) <- momentNames
    return(momentMatrix)
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

#' Print results
#'
#' This function uses the print method on the ivmte return list.
#'
#' @param x an object returned from '\code{ivmte}'.
#' @param ... additional arguments.
#' @return basic set of results.
#' @export
print.ivmte <- function(x, ...) {
    stopifnot(inherits(x, "ivmte"))
    if (!is.null(x$bounds)) {
        cat("\n")
        ## Return bounds, audit cout, and minumum criterion
        cat(sprintf("Bounds on the target parameter: [%s, %s]\n",
                    fmtResult(x$bounds[1]),
                    fmtResult(x$bounds[2])))
        if (!is.null(x$audit.grid$violations)) {
            cat(sprintf("Audit reached audit.max (%s)\n",
                        x$audit.count))
        } else {
            if (x$audit.count == 1) rs <- "round"
            if (x$audit.count > 1) rs <- "rounds"
            cat(sprintf("Audit terminated successfully after %s",
                        x$audit.count), rs, "\n")
        }
        cat(sprintf("Minimum criterion: %s \n", x$audit.criterion))
        ## Return LP solver used
        cat(sprintf("LP solver: %s\n", x$lpsolver))
        if (x$lpsolver == "lp_solve ('lpSolveAPI')") {
            warning(
                gsub("\\s+", " ",
                     "The R package 'lpSolveAPI' interfaces with 'lp_solve',
                      which is outdated and potentially unreliable.  It is
                      recommended to use commercial solvers
                      Gurobi (lpsolver = 'gurobi')
                      or CPLEX (lpsolver = 'cplexAPI') instead.
                      Free academic licenses can be obtained for these
                      commercial solvers."),
                "\n", call. = FALSE)
        }
    }
    if (!is.null(x$pointestimate)) {
        cat("\n")
        ## Return point estimate
        cat(sprintf("Point estimate of the target parameter: %s\n",
                    fmtResult(x$pointestimate)))
    }
    cat("\n")
}

#' Summarize results
#'
#' This function uses the summary method on the ivmte return list.
#'
#' @param object an object returned from '\code{ivmte}'.
#' @param ... additional arguments.
#' @return summarized results.
#' @export
summary.ivmte <- function(object, ...) {
    stopifnot(inherits(object, "ivmte"))
    ## Summary for the partially identified case
    if (!is.null(object$bounds)) {
        cat("\n")
        ## Return bounds, audit count, and minumum criterion
        cat(sprintf("Bounds on the target parameter: [%s, %s]\n",
                    fmtResult(object$bounds[1]),
                    fmtResult(object$bounds[2])))
        if (!is.null(object$audit.grid$violations)) {
            cat(sprintf("Audit reached audit.max (%s)\n",
                        object$audit.count))
        } else {
            if (object$audit.count == 1) rs <- "round"
            if (object$audit.count > 1) rs <- "rounds"
            cat(sprintf("Audit terminated successfully after %s",
                        object$audit.count), rs, "\n")
        }
        cat(sprintf("Minimum criterion: %s \n", object$audit.criterion))
        ## Return LP solver used
        cat(sprintf("LP solver: %s\n", object$lpsolver))
        if (object$lpsolver == "lp_solve ('lpSolveAPI')") {
            warning(
                gsub("\\s+", " ",
                     "The R package 'lpSolveAPI' interfaces with 'lp_solve',
                      which is outdated and potentially unreliable.  It is
                      recommended to use commercial solvers
                      Gurobi (lpsolver = 'gurobi')
                      or CPLEX (lpsolver = 'cplexAPI') instead.
                      Free academic licenses can be obtained for these
                      commercial solvers."),
                "\n", call. = FALSE)
        }
        if (!is.null(object$bootstraps)) {
            ## Return confidence intervals and p-values
            levels <- as.numeric(rownames(object$bounds.ci[[1]]))
            ## ciTypes <- names(object$bounds.ci)
            ciTypes <- "backward"
            ciN <- 1
            for (i in ciTypes) {
                ## Confidence intervals
                cat("\nBootstrapped confidence intervals (",
                    i, "):\n", sep = "")
                for (j in 1:length(levels)) {
                    cistr <- paste0("[",
                                    fmtResult(object$bounds.ci[[i]][j, 1]),
                                    ", ",
                                    fmtResult(object$bounds.ci[[i]][j, 2]),
                                    "]")
                    cat("    ",
                        levels[j] * 100,
                        "%: ",
                        cistr, "\n", sep = "")
                }
                ## p-values
                cat("p-value: ",
                    fmtResult(object$pvalue[ciN]), "\n", sep = "")
                ciN <- ciN + 1
            }
            ## Return specification test
            if (!is.null(object$specification.pvalue)) {
                cat("Bootstrapped specification test p-value: ",
                    fmtResult(object$specification.pvalue), "\n", sep = "")
            }
            ## Return bootstrap counts
            cat(sprintf("Number of bootstraps: %s",
                        object$bootstraps))
            if (object$bootstraps.failed > 0) {
                cat(sprintf(" (%s failed and redrawn)\n",
                            object$bootstraps.failed))
            } else {
                cat("\n")
            }
        }
    }
    ## Summary for the point identified case
    if (!is.null(object$pointestimate)) {
        cat("\n")
        ## Return bounds, audit cout, and minumum criterion
        cat(sprintf("Point estimate of the target parameter: %s\n",
                    fmtResult(object$pointestimate)))
        if (!is.null(object$bootstraps)) {
            ## Return confidence intervals and p-values
            levels <- as.numeric(rownames(object$pointestimate.ci[[1]]))
            ciTypes <- c("nonparametric", "parametric")
            for (i in 1:1) { ## Note: 1:1 is deliberate, only want to
                             ## present nonparametric CI
                if (i == 1) ciType <- "nonparametric"
                if (i == 2) ciType <- "normal quantiles"
                ## Confidence intervals
                cat("\nBootstrapped confidence intervals (",
                    ciType, "):\n", sep = "")
                for (j in 1:length(levels)) {
                    cistr <-
                        paste0("[",
                               fmtResult(object$pointestimate.ci[[i]][j, 1]),
                               ", ",
                               fmtResult(object$pointestimate.ci[[i]][j, 2]),
                               "]")
                    cat("    ",
                        levels[j] * 100,
                        "%: ",
                        cistr, "\n", sep = "")
                }
                ## p-values
                cat("p-value: ",
                    fmtResult(object$pvalue[i]), "\n", sep = "")
            }
            ## Return specification test
            if (!is.null(object$jtest)) {
                cat("Bootstrapped J-test p-value: ",
                    fmtResult(object$jtest[2]), "\n", sep = "")
            }
            ## Return bootstrap counts
            cat(sprintf("Number of bootstraps: %s",
                        object$bootstraps))
            if (object$bootstraps.failed > 0) {
                cat(sprintf(" (%s failed and redrawn)\n",
                            object$bootstraps.failed))
            } else {
                cat("\n")
            }
        }
    }
    cat("\n")
}
