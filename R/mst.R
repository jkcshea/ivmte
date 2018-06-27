#' Estimation procedure from Mogstad, Torgovitsky (2017)
#'
#' This function estimates bounds on treatment effect parameters,
#' following the procedure described in Mogstad, Torgotvitsky
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
#'     correpsonding to each IV-like estimand. See
#'     \code{\link[mst]{list.mst}} on how to input the argument.
#' @param components a list of vectors of the terms/components from
#'     the regressions specifications we want to include in the set of
#'     IV-like estimands. To select the intercept term, include in the
#'     vector of variable names, `intercept'. See
#'     \code{\link[mst]{list.mst}} on how to input the argument.
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
#' @param uname variable name for unobservale used in declaring MTRs.
#' @param target target parameter to be estimated. Currently function
#'     allows for ATE ("ate"), ATT ("att"), ATU ("atu"), LATE
#'     ("late"), and generalized LATE ("genlate").
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
#' @param audit.Nx number of points in the interval [0, 1],
#'     corresponding to the normalized value of the unobservable term,
#'     to audit in each iteration of the audit procedure.
#' @param audit.max maximum number of iterations in the audit
#'     procedure.
#' @param audit.tol tolerance for determining when to end the audit
#'     procedure.
#' @param m1.ub numeric value for upper bound on MTR for treated
#'     group.
#' @param m0.ub numeric value for upper bound on MTR for control
#'     group.
#' @param m1.lb numeric value for lower bound on MTR for treated
#'     group.
#' @param m0.lb numeric value for lower bound on MTR for control
#'     group.
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
#' @return Returns a list of results from throughout the estimation
#'     procedure. This includes all IV-like estimands; the propensity
#'     score model; bounds on the treatment effect; the estimated
#'     expectations of each term in the MTRs; the components and
#'     results of the LP problem.
#'
#' @examples
#' ivlikespecs <- c(ey ~ d | z,
#'                ey ~ d | factor(z),
#'                ey ~ d,
#'                ey ~ d | factor(z))
#' jvec <- lists.mst(d, d, d, d)
#' svec <- lists.mst(, , , z %in% c(2, 4))
#'
#' mst(ivlike = ivlikespecs,
#'     data = dtm,
#'     components = jvec,
#'     propensity = pz,
#'     subset = svec,
#'     m0 = ~  u + u^2,
#'     m1 = ~  u + u^2,
#'     uname = u,
#'     target = "att",
#'     m0.dec = TRUE,
#'     m1.dec = TRUE)
#'
#'
#' mst(ivlike = y ~ d + x1 + x2 | x1 + x2 + z1 + z2, dt,
#'     components = d,
#'     propensity = d ~ z1 + z2 + x1 + x2,
#'     link = "logit",
#'     m0 = ~ x1 + x2,
#'     m1 = ~ x1  ,
#'     target = "late",
#'     Z = c(z1, z2),
#'     from = c(0,1),
#'     to= c(1,2))
#'
#' @export
mst <- function(ivlike, data, subset, components, propensity,
                link, treat, m0, m1, uname = u, target, late.Z,
                late.from, late.to, late.X, eval.X, genlate.lb, genlate.ub,
                obseq.tol = 1, grid.Nu = 20, grid.Nx = 20,
                audit.Nx = 2, audit.Nu = 3, 
                audit.max = 5, audit.tol = 1e-08,
                m1.ub, m0.ub, m1.lb, m0.lb, mte.ub,
                mte.lb, m0.dec, m0.inc, m1.dec, m1.inc, mte.dec,
                mte.inc) {

    ## Match call arguments
    call     <- match.call(expand.dots = FALSE)

    ## FIX: at some point, you may want to include "weights"

    ##---------------------------
    ## 0.a Check format of `formula', `subset', and `component' inputs
    ##---------------------------
    
    if (class_list(ivlike)) {
        
        ## Convert formula, components, and subset inputs into lists
        length_formula <- length(ivlike)

        if (hasArg(components)) {
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
                         component vectors. Specifications without coresponding
                         component vectors will include all covariates when
                         constructing the S-set."),
                    call. = FALSE)
            components[(length(components) + 1) : length(ivlike)] <- ""
        }
        
        if (length_formula < length_components) {
            warning(gsub("\\s+", " ",
                         "List of components not the same length of list of
                         IV-like specifications: more component vectors than
                         specifications. Component vectors without coresponding
                         specifications will be dropped."),
                    call. = FALSE)
            components <- components[1 : length(ivlike)]
        }

        ## Check the subset input---of the three lists that are input,
        ## only this can be omitted by the user, in which case no
        ## subsetting is used
        if (hasArg(subset)) {
            if (!class(subset) == "list") subset <- list(subset)
            if(length(subset) > length_formula) {
                warning(gsub("\\s+", " ",
                             "List of subset conditions not the same length
                              of list IV-like specifications: more
                              specifications than subsetting conditions.
                              Specifications without corresponding subset
                              conditions will include all observations."),
                        call. = FALSE)
                subset[length(subset) + 1 : length(ivlike)] <- ""
            }
            if(length(subset) < length_formula) {
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
    ## 0.b Check monotonicity conditions
    ##---------------------------

    if (hasArg(m0.inc) & hasArg(m0.dec)) {
        if (m0.inc == TRUE & m0.dec == TRUE) {
            stop(gsub("\\s+", " ",
                      "Cannot have m0 be monotone increasing and monotone
                      decreasing."))
        }
    }

    if (hasArg(m1.inc) & hasArg(m1.dec)) {
        if (m1.inc == TRUE & m1.dec == TRUE) {
            stop(gsub("\\s+", " ",
                      "Cannot have m1 be monotone increasing and monotone
                      decreasing."))
        }
    }

    if (hasArg(mte.inc) & hasArg(mte.dec)) {
        if (mte.inc == TRUE & mte.dec == TRUE) {
            stop(gsub("\\s+", " ",
                      "Cannot have MTE be monotone increasing and monotone
                      decreasing."))
        }
    }

    ##---------------------------
    ## 0.c Check numeric arguments and case completion
    ##---------------------------

    if (hasArg(treat)) {
        if (! deparse(substitute(treat)) %in% colnames(data)) {
            stop("Declared treatment indicator not found in data")
        }
    }

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
                      "'genlate.lb' and 'genlate.ub' must be between 0 and 1."))
        }
        if (genlate.lb >= genlate.ub) {
            stop("'genlate.lb' must be less than 'genlate.ub'.")
        }
    }
    if (target != "genlate" & (hasArg("genlate.lb") | hasArg("genlate.ub"))) {
        warning(gsub("\\s+", " ",
                     "Unless target parameter is 'genlate', 'genlate.lb' and
                     'genlate.ub' arguments will not be used."))
    }
    if (target != "late" & (hasArg("eval.X") | hasArg("late.X"))) {
        warning(gsub("\\s+", " ",
                     "Unless target parameter is 'late', 'eval.X' and 'late.X'
                     arguments will not be used."))
    }
    if (hasArg(link)) {
        if (! link %in% c("linear", "logit", "probit")) {
            stop(gsub("\\s+", " ",
                      "Specified link is not recognized. Choose from 'linear',
                      'logit', or 'probit'."))
        }
    }
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

    ## FIX: you need to remove the spline terms, but you also need to
    ## keep track of the variables that interact with the
    ## splines. This should not be hard since the list you generate
    ## includes the variables that interact with the spline.
    
    ## Restrict data used for all parts of procedure to be the same.
    ## Collect list of all terms used in formula

    vars_y          <- NULL
    vars_formulas_x <- c()
    vars_formulas_z <- c()
    vars_subsets    <- c()
    vars_mtr        <- c()
    vars_propensity <- c()

    terms_formulas_x <- c()
    terms_formulas_z <- c()
    terms_mtr        <- c()
    
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

        terms_mtr <- c(attr(terms(splinesobj[[1]]$formula), "term.labels"),
                       attr(terms(splinesobj[[2]]$formula), "term.labels"))

        if (!is.null(splinesobj[[1]]$splineslist)) {
            sf0 <- as.formula(paste("~",
                                    paste(unlist(splinesobj[[1]]$splineslist),
                                          collapse = " + ")))
            vars_mtr <- c(vars_mtr, all.vars(sf0))
            terms_mtr <- c(terms_mtr, attr(terms(sf0), "term.labels"))
        }
      
        if (!is.null(splinesobj[[2]]$splineslist)) {
            sf1 <- as.formula(paste("~",
                                    paste(unlist(splinesobj[[2]]$splineslist),
                                          collapse = " + ")))
            vars_mtr <- c(vars_mtr, all.vars(sf1))
            terms_mtr <- c(terms_mtr, attr(terms(sf1), "term.labels"))
        }
        
    } else {
        stop("m0 and m1 must be one-sided formulas.")
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
                          "Propensity score argument is interpretted as a
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
                             treatment variable."))
                treat <- all.vars(ivlike)[2]
                vars_propensity <- c(vars_propensity,
                                      treat)
            } else if (is.list(ivlike)) {
                warning(gsub("\\s+", " ",
                             "First independent variable of first IV regression
                             is selected as the treatment variable."))
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
                              unlist(terms_mtr))

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
                    terms_prxopensity)
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
                         "Specifications without coresponding
                         component vectors will include all covariates when
                         constructing the S-set."),
                    call. = FALSE)
        }
        components[compMissing] <- comp_filler[compMissing]
    } else {
        components <- comp_filler

    }  
   
    ## Keep only complete cases

    varError <- allvars[! allvars %in% colnames(data)]
    if (length(varError) > 0) {
        varError <- paste0("The following variables are not contained
                          in the data set: ", varError, ".")
        stop(gsub("\\s+", " ", varError), call. = FALSE)
    }
    
    data  <- data[(complete.cases(data[, allvars])), ]
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
    
    ## Generate target S- weights
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
        stop("Unrecognized target parameter (or custom weights not ready)")
    }

    ## Integrate m0 and m1 functions
    m0call <- modcall(call,
                      newcall = polyparse.mst,
                      keepargs = c("uname"),
                      newargs = list(formula = m0,
                                     data = quote(cdata)))

    m1call <- modcall(call,
                      newcall = polyparse.mst,
                      keepargs = c("uname"),
                      newargs = list(formula = m1,
                                     data = quote(cdata)))

    pm0 <- eval(as.call(m0call))
    pm1 <- eval(as.call(m1call))

    ## Estimate gamma-star
    gstar0 <- gengamma.mst(pm0, w0$lb, w0$ub, w0$mp)
    gstar1 <- gengamma.mst(pm1, w1$lb, w1$ub, w1$mp)
  
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
   
    gstar0 <- c(gstar0, gstarSpline0)
    gstar1 <- c(gstar1, gstarSpline1)

    ##---------------------------
    ## 4. Generate moments/gamma terms for IV-like estimands
    ##---------------------------

    message("Generating IV-like moments...")
    
    sset  <- list() ## Contains all IV-like estimates and their
                    ## coresponding moments/gammas
    scount <- 1     ## counter for S-set constraints

    ## Construct `sset' object when a single IV-like specification is
    ## provided
    if (class_formula(ivlike)) {

        ## Obtain coefficient estimates and S-weights
        scall <- modcall(call,
                       newcall = sweights.mst,
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
            sest  <- sweights.mst(formula = sformula,
                                  data = sdata,
                                  component = scomponent,
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
                                         terms_mtr = quote(terms_mtr),
                                         sset = quote(sset),
                                         gstar0 = quote(gstar0),
                                         gstar1 = quote(gstar1)))

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
                lpresult =  audit$lpresult))
}


#' Generating LP moments for IV-like estimands
#'
#' This function takes in the IV estimate and its IV-like
#' specification, and generates a list containing the corresponding
#' point estimate, and the corresponding moments (gammas) that will
#' enter into the contraint matrix of the LP problem.
#'
#' @param sset A list, which is modified and returned as the output.
#' @param sest A list containing the point estimates and S-weights
#'     corresponding to a particular IV-like estimand.
#' @param pmodobj A vector of propensity scores.
#' @param pm0 A list of the monomials in the MTR for d = 0.
#' @param pm1 A list of the monomials in the MTR for d = 1.
#' @param ncomponents The number of components from the IV regression
#'     we want to include in the S-set.
#' @param scount A counter for the number of elements in the S-set.
#' @param subset_index An index for the subset of the data the IV
#'     regression is restricted to.
#' @return A list containing the point estimate for the IV regression,
#'     and the expectation of each monomoial term in the MTR.
gensset.mst <- function(data, sset, sest, splinesobj, pmodobj, pm0, pm1,
                        ncomponents, scount, subset_index) {

    for (j in 1:ncomponents) {
        message(paste0("    Moment ", scount, "..."))
        gs0 <- gengamma.mst(monomials = pm0,
                            lb = pmodobj,
                            ub = 1,
                            multiplier = sest$sw0[, j],
                            subset = subset_index)

        gs1 <- gengamma.mst(monomials = pm1,
                            lb = 0,
                            ub = pmodobj,
                            multiplier = sest$sw1[, j],
                            subset = subset_index)
        
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

