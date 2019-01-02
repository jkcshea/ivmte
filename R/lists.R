#' Listing subsets and components
#'
#' This function allows the user to declare a list of variable names
#' in non-character form and subsetting conditions. This is used to
#' ensure clean entry of arguments into the \code{components} and
#' \code{subset} arguments of the function.
#'
#' @param ... subset conditions or variable names
#' @return list.
#'
#' @examples
#' components <- l(d, x1, intercept)
#' subsets <- l(, z %in% c(2, 4))
#'
#' @export
l <- function(...) {
    lists <- as.list(substitute(list(...)))
    return(lists[-1])
}
