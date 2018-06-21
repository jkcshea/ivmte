#' Auxiliary function: \code{which} for lists
#' 
#' Auxiliary function that makes it possible to use \code{which} with
#' a list.
#' @param vector the vector for which we want to check the entries of
#' @param obj the value for which we want the vector to match on.
#' @return a vector of positions where the elements in \code{vector}
#'     are equal to \code{obj}.
whichforlist <- function(vector, obj) {
    which(vector == obj)
}

#' Auxiliary function: extracting elements from strings
#' 
#' This auxiliary function extracts the (string) element in the
#' \code{position} argument of the \code{vector} argument.
#' @param vector the vector from which we want to extract the
#'     elements.
#' @param position the position in \code{vector} to extract.
#' @param truncation the number of characters from the front of the
#'     element being extracted that should be dropped.
#' @return A chracter/string.
vecextract <- function(vector, position, truncation = 0) {
    elem <- vector[position]
    elem <- substr(elem, truncation + 1, nchar(elem))
    return(elem)
}

#' Generating monomials
#' 
#' This function takes in a first vector of coefficients, and a second
#' vector declaring which univariate polynomial basis corresponds to
#' each element of the first vector. Then it generates a list of
#' monomials corresponding to the polynomial.
#' @param vector vector of polynomial coefficients.
#' @param basis vector of exponents corresponding to each element of
#'     \code{vector}.
#' @param zero logical, if \code{FALSE} then \code{vector} does not
#'     include an element for the constant term. The vector
#'     \code{basis} will need to be adjusted to account for this in
#'     order to generate the correct polynomial and monomials.
#' @return A list of monomials, in the form of the \code{polynom}
#'     package.
genmono <- function(vector, basis, zero = FALSE) {
    if (!zero) basis <- basis + 1
    monolist  <- mapply(genej, pos = basis, length = basis)
    polyinput <- mapply("*", monolist, vector)
    poly <- lapply(polyinput, polynom::polynomial)    
    return(poly)
}


#' Evaluating polynomials
#' 
#' This function allows one to evaluate a list of polynomials at
#' various points (the points are allowed to differ across
#' polynomials). This function will instead be used to evaluate a list
#' of monomials.
#' @param polynomials a list of polynomials.
#' @param points a list/vector of points at which we want to evaluate
#'     each polynomial.
#' @return A matrix of values, corresponding to the polynomials
#'     specified in \code{polynomials} evaluated at the points
#'     specified in \code{points}.
polylisteval <- function(polynomials, points) {
  
    if (length(polynomials) != length(points)) {
        stop(gsub("\\s+", " ",
                  "List of polynomials to evaluate, and list of points to
                  evaluate each polynomial, are not equal."))
    }
    output <- mapply(predict, polynomials, points)
    return(output)
}

#' Parsing marginal treatment response formulas
#' 
#' This function takes in an MTR formula, and then parses the formula
#' such that it becomes a polynomial in the unobservable u. It then
#' breaks these polynomials into monomials, and then integrates each
#' of them with respect to u. Each integral corresponds to E[md | D,
#' X, Z].
#' @param formula the MTR.
#' @param data \code{data.frame} for which we obtain E[md | D, X, Z] for each
#'     observation.
#' @param uname variable name for unobservable used in declaring the
#'     MTR.
#' @return A list (of lists) of monomials corresponding to the
#'     original MTR (for each observation); a list (of lists) of the
#'     integrated monomials; a vector for the degree of each of the
#'     original monomials in the MTR; and a vector for the names of
#'     each variable entering into the MTR (note \code{x^2 + x} has
#'     only one term, \code{x}).
polyparse.mst <- function(formula, data, uname = u) {

    ## FIX: you need to include a component that will separate out the
    ## splines. If you allow for all the arguments to be put into the
    ## fomrula, that seems like you'd get some hideous outcome. 
    
    ## update formula parsing
    formula <- Formula::as.Formula(formula)
    uname   <- deparse(substitute(uname))
    
    ## Include redundant variable u, so monomials in m0, m1
    ## specifications correspond polynomial coefficients on u
    ## monomials
    data[[uname]] <- 1
    dmat <- design.mst(formula, data)$X 
    
    ## Save original terms for reference as later
    oterms <- attr(terms(formula), "term.labels")

    ## Separate and parse the original terms
    nterms <- oterms
    for (sep in c("I\\(", "\\)", "\\*", ":", "  ", "  ")) {
        nterms <- lapply(nterms, gsub, pattern = sep, replacement = " ")
    }
    nterms <- strsplit(trimws(unlist(nterms)), " ")

    ## Find monomials of degree 1 ('degree' is with respect to u)
    u_pos <- lapply(nterms, whichforlist, obj = uname)
    u_pos <- which(u_pos > 0)
    
    ## Find monomials of degree exceeding 1, and determine their degree
    trunc_nterms <- lapply(nterms, substr, 0, nchar(uname) + 1)
    uexp_pos     <- lapply(trunc_nterms, whichforlist, obj = paste0(uname, "^"))
    uexp_pos     <- which(uexp_pos > 0)
    uexp_subpos  <- lapply(trunc_nterms[uexp_pos],
                           whichforlist, obj =  paste0(uname, "^"))
    deggtr2      <- FALSE
    if (length(unlist(uexp_subpos)) > 0) deggtr2 <- TRUE    
    
    ## Create a matrix stating the degree for each monomial with degree >= 1
    if(length(u_pos) == 0){
        exptab1 <- NULL
    } else {
        exptab1 <- cbind(u_pos, 1)
    }
    if (deggtr2) {
        uexp <- as.numeric(mapply(vecextract,
                                  nterms[uexp_pos],
                                  position = uexp_subpos,
                                  truncation = nchar(uname) + 1))
        exptab2 <- cbind(uexp_pos, uexp)
        exptab  <- rbind(exptab1, exptab2)
    } else {
        exptab <- exptab1
    }

    if(!is.null(exptab)) {
        colnames(exptab) <- c("term", "degree")
    }
    
    ## Determine which terms do not involve u
    nonuterms    <- unlist(oterms[!seq(1, length(oterms)) %in% exptab[, 1]])
    nonutermspos <- which(oterms %in% nonuterms)
    if (length(nonuterms) > 0) {
        exptab0 <- cbind(nonutermspos, replicate(length(nonutermspos), 0))
    } else {
        exptab0  <- NULL
    }
    exptab <- rbind(exptab0, exptab)
    exptab <- exptab[order(exptab[,1]),]
    if (is.matrix(exptab)) {
        exporder <- exptab[, 2]
        colnames(exptab) <- c("term", "degree")
    } else {
        exporder <- exptab[2]
        names(exptab) <- c("term", "degree")
    }
    names(exporder) <- NULL
    
    ## generate matrix with monomial coefficients
    if ("(Intercept)" %in% colnames(dmat)) {
        exporder <- c(0, exporder)
        oterms   <- c("(Intercept)", oterms)
    }
    polymat <- dmat[, oterms]
    
    ## prepare monomials and their integrals
    monomial_list <- lapply(split(polymat, seq(nrow(polymat))),
                            genmono,
                            basis = exporder)

    integral_list <- lapply(monomial_list,
                            lapply,
                            polynom::integral)

    names(monomial_list) <- rownames(data)
    names(integral_list) <- rownames(data)
    
    return(list(mlist = monomial_list,
                ilist = integral_list,
                exporder = exporder,
                terms = oterms))  

    ## Note: to implement numerical integration, you should be able to
    ## take these polynomials you have parsed, multiply them by the
    ## IV-like specifications, and then apply numerical integration.
}

#' Estimating expectations of terms in the MTR (gamma objects)
#' 
#' This function generates the gamma objects defined in the paper,
#' i.e. each additive term in E[md], where md is a MTR.
#'
#' @param monomials object containing list of list of monomials. Each
#'     element of the outer list represents an observation in the data
#'     set, each element in the inner list is a monomial from the
#'     MTR. The variable is the unobservable u, and the coefficient is
#'     the evaluation of any interactions with u.
#' @param lb vector of lower bounds for the interval of
#'     integration. Each element corresponds to an observation.
#' @param ub vector of upper bounds for the interval of
#'     integration. Each element corresponds to an observation.
#' @param multiplier a vector of the weights that enter into the
#'     interal. Each element corresponds to an observation.
#' @param subset Subset condtion used to select observations with
#'     which to estimate gamma.
#' @param means logical, if TRUE then function returns the terms of
#'     E[md]. If FALSE, then function instead returns each term of
#'     E[md | D, X, Z]. This is useful for testing the code,
#'     i.e. obtaining population estimates.
#' @return If \code{means = TRUE}, then the function returns a vector
#'     of the additive terms in Gamma (i.e. the expectation is over D,
#'     X, Z, and u). If \code{means = FALSE}, then the function
#'     returns a matrix, where each row corresponds to an observation,
#'     and each column corresponds to an addiive term in E[md | D, X,
#'     Z] (i.e. only the integral with respect to u is performed).
#'
#' @export 
gengamma.mst <- function(monomials, lb, ub, multiplier = 1,
                         subset = NULL, means = TRUE) {

    ## FIX: You need to allow for a spline components here.
    
    exporder  <- monomials$exporder
    integrals <- monomials$ilist
    
    if(!is.null(subset)) integrals <- integrals[subset]
    nmono     <- length(exporder)

    ## Determine bounds of integrals (i.e. include weights)
    if (length(ub) == 1) ub <- replicate(length(integrals), ub)
    if (length(lb) == 1) lb <- replicate(length(integrals), lb)
    
    ub <- split(replicate(nmono, ub), seq(length(ub)))
    lb <- split(replicate(nmono, lb), seq(length(lb)))
   
    monoeval <- t(mapply(polylisteval, integrals, ub)) -
        t(mapply(polylisteval, integrals, lb))
    
    if (means) {
        gstar <- colMeans(monoeval * multiplier)
        names(gstar) <- monomials$terms

        return(gstar)
    } else {
        return(monoeval * multiplier)
    }
}

#' Separating splines from MTR formulas
#'
#' This function separates out the function call "uSpline()"
#' potentially embedded in the MTR formulas from the rest of the
#' fomrula. The terms involving splines are treated separately from
#' the terms that do not involve splines when creating the gamma
#' moments.
#' @param formula the formula that is to be parsed.
#' @return a list containing two objects. One object is \code{formula}
#'     but with the spline components removed. The second object is a
#'     list. The name of each element is the "uSpline()" command, and
#'     the elements are a vector of the names of covariates that were
#'     interacted with the "uSpline()" command.
removeSplines <- function(formula) {

    fterms <- attr(terms(formula), "term.labels")
    whichspline <- sapply(fterms,
                          function(y) grepl(x = y, pattern = "uSpline\\("))
    if (max(whichspline) == 1) {
        ftobj <- terms(formula)
        splinepos <- which(whichspline == TRUE)
        nosplines <- drop.terms(ftobj, splinepos)
        nosplines <- Formula::as.Formula(nosplines)
        
        splineterms <- fterms[whichspline]
        splinelist <- list()
        for (splineobj in splineterms) {
            
            splinepos <- regexpr("uSpline\\(", splineobj)
            degreepos <- regexpr("degree = ", splineobj)
            knotspos  <- regexpr("knots = ", splineobj)
            knotslpos <- regexpr("knots = c\\(", splineobj)
            interpos  <- regexpr("intercept = ", splineobj)

            substr(splineobj, splinepos, nchar(splineobj))

            ## Check if vectors or sequences are declared to adjust parsing
            firstopen  <- regexpr("\\(",
                                  substr(splineobj,
                                         splinepos + 8,
                                         nchar(splineobj)))
            firstclose <- regexpr("\\)",
                                  substr(splineobj,
                                         splinepos,
                                         nchar(splineobj)))
            secondclose <- regexpr("\\)",
                                   substr(splineobj,
                                          splinepos + 8 + firstclose,
                                          nchar(splineobj)))

            if ((firstopen < firstclose) & (firstopen != - 1)) {
                splinecmd <- substr(splineobj,
                                    splinepos,
                                    splinepos + 8 + firstclose + secondclose)
            } else {
                splinecmd <- substr(splineobj,
                                    splinepos,
                                    splinepos + firstclose - 1)
            }

            ## Separate uSpline command from terms interacting with the spline
            splinecmdstr <- gsub("\\)", "\\\\)",
                                 gsub("\\(", "\\\\(", splinecmd))   

            interobj <- gsub(paste0(":", splinecmdstr), "", splineobj)
            interobj <- gsub(paste0(splinecmdstr, ":"), "", interobj)
            interobj <- gsub("::", ":", interobj)
            
            if (interobj == splinecmd) interobj <- "1"
            
            if (splinecmd %in% names(splinelist)) {
                splinelist[[splinecmd]] <- c(splinelist[[splinecmd]], interobj)
            } else {
                splinelist[[splinecmd]] <- c(interobj)
            }
        }
        
    } else {
        nosplines <- formula
        splinelist <- NULL
    }
    
    return(list(formula = nosplines,
                splinelist = splinelist))
}

#' Integrated splines
#'
#' This function integrates out splines that the user specifies when
#' declaring the MTRs. This is to be used when generating the gamma
#' moments.
#' @param x the points to evluate the interal of the the splines.
#' @param knots the knots of the spline.
#' @param degree the degree of the spline; default is set to 0
#'     (constant splines).
#' @param intercept boolean, set to TRUE if intercept term is to be
#'     included (i.e. an additional basis such that the sum of the
#'     splines at every point in \code{x} is equal to 1).
#' @return a matrix, the values of the integrated splines. Each row
#'     corresponds to a value of \code{x}; each column corresponds to
#'     a basis defined by the degrees and knots.
#'
#' @examples
#' Since the splines are declared as part of the MTR, you will need to
#' have parsed out the spline command. Thus, this command will be
#' called via eval(parse(text = .)). In the examples below, the
#' commands are parsed from the object "splinelist" generated by
#' \code{\link[MST]{removeSplines}}. The names of the elements in the
#' list are the spline commands, and the elements themselves are the
#' terms that interact with the splines.
#'
#' eval(parse(text = gsub("uSpline\\(",
#'                        "uSplineInt(x = x, ",
#'                        names(splinelist)[1])))
#'
#' eval(parse(text = gsub("uSpline\\(",
#'                        "uSplineInt(x = x, ",
#'                         names(splinelist)[2])))
uSplineInt <- function(x, knots, degree = 0, intercept = TRUE) {
    ibs(x = x,
        knots = knots,
        degree = degree,
        intercept = intercept,
        Boundary.knots = c(0, 1))
}

#' Spline basis function
#'
#' This function evaluates the splines that the user specifies when
#' declaring the MTRs. This is to be used for auditing, namely when
#' checking the boundedness and monotonicity conditions.
#' @param x the points to evluate the interal of the the splines.
#' @param knots the knots of the spline.
#' @param degree the degree of the spline; default is set to 0
#'     (constant splines).
#' @param intercept boolean, set to TRUE if intercept term is to be
#'     included (i.e. an additional basis such that the sum of the
#'     splines at every point in \code{x} is equal to 1).
#' @return a matrix, the values of the integrated splines. Each row
#'     corresponds to a value of \code{x}; each column corresponds to
#'     a basis defined by the degrees and knots.
#'
#' @examples
#'
#' Since the splines are declared as part of the MTR, you will need to
#' have parsed out the spline command. Thus, this command will be
#' called via eval(parse(text = .)). In the examples below, the
#' commands are parsed from the object "splinelist" generated by
#' \code{\link[MST]{removeSplines}}. The names of the elements in the
#' list are the spline commands, and the elements themselves are the
#' terms that interact with the splines.
#' 
#' eval(parse(text = gsub("uSpline\\(",
#'                        "uSplineBasis(x = x, ",
#'                         names(splinelist)[1])))
#'
#' eval(parse(text = gsub("uSpline\\(",
#'                        "uSplineBasis(x = x, ",
#'                        names(splinelist)[2])))
uSplineBasis <- function(x, knots, degree = 0, intercept = TRUE) {
    bSpline(x = x,
            knots = knots,
            degree = degree,
            intercept = intercept,
            Boundary.knots = c(0, 1))
}
