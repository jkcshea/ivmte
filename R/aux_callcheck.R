#' Auxiliary function: remove extraneous spaces
#'
#' Auxiliary function to remove extraneous spaces from strings.
#' @param string the string object to be cleaned.
#' @return a string
subsetclean <- function(string) {
    string <- gsub(", ", ",", string)
    string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE)
    return(unlist(strsplit(string, split = " ")))
}

#' Auxiliary function: check if string is command
#'
#' Auxiliary function to check if a string is in fact a command, but
#' in string form.
#' @param string the string object to be checked.
#' @return boolean expression.
isfunctionstring <- function(string) {
    regexpr("\\(", string) > 0 & substr(string, nchar(string), nchar(string)) == ")"
}

#' Auxiliary function: extract arguments from function in string form
#'
#' Auxiliary function to extract arguments from a function that is in
#' string form.
#' @param string the function in string form.
#' @return string of arguments.
argstring <- function(string) {
    pos <- regexpr("\\(", string)
    args <- substr(string, pos + 1, nchar(string) - 1)
    return(args)
}

#' Auxiliary function: test if object is a formula
#'
#' Auxiliary function to test if an object is a formula. Warnings are
#' suppressed.
#' @param obj the object to be checked.
#' @return boolean expression.
class_formula <- function(obj) {
    suppressWarnings(try(class(obj), silent = TRUE) == "formula")
}

#' Auxiliary function: test if object is a list
#'
#' Auxiliary function to test if an object is a list. Warnings are
#' suppressed.
#' @param obj the object to be checked.
#' @return boolean expression.
class_list <- function(obj) {
    suppressWarnings(try(class(obj), silent = TRUE) == "list")
}
    
#' Auxiliary function: extract X and Z covariates from a formula
#'
#' Auxiliary function that takes in a formula, and extracts the
#' variable names of either the covariates or instruments.
#' @param fm the formula.
#' @param inst boolean expression, set to TRUE if the instrument names
#'     are to be extracted. Otherwise, the covariate names are
#'     extracted.
get_xz <- function(fm, inst = FALSE) {
    fm <- Formula::as.Formula(fm)
    if (length(fm)[2] == 1) {
        x <- all.vars(fm)[-1]
        z <- NULL
    }
    if (length(fm)[2] == 2) {
        x <- all.vars(formula(fm, rhs = 1))[-1]
        z <- all.vars(formula(fm, rhs = 2))[-1]
    }
    if (inst == TRUE) {
        return(z)
    } else {
        return(x)
    }
}
