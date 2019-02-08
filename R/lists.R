#' Listing subsets and components
#'
#' This function allows the user to declare a list of variable names
#' in non-character form and subsetting conditions. This is used to
#' ensure clean entry of arguments into the \code{components} and
#' \code{subset} arguments of the function. When selecting components
#' to include in the S set, selecting the intercept term and factor
#' variables requires special treatment. To select the intercept term,
#' include in the vector of variable names, `intercept'. If the the
#' factorized counterpart of a variable \code{x = 1, 2, 3} is included
#' in the IV-like specifications via \code{factor(x)}, the user can
#' select the coefficients for specific factors by declaring the
#' components \code{factor(x)-1, factor(x)-2, factor(x)-3}.
#'
#' @param ... subset conditions or variable names
#' @return list.
#'
#' @examples
#' components <- l(d, x1, intercept, factor(x)-2)
#' subsets <- l(, z %in% c(2, 4))
#'
#' @export
l <- function(...) {
    lists <- as.list(substitute(list(...)))
    return(lists[-1])
}
