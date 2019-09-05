#' AngristEvans Data
#'
#' @format A data frame with 209,133 rows and 5 columns.
#' \describe{
#'   \item{worked}{indicator for whether worked in the previous year}
#'   \item{hours}{weekly hours worked in the previous year}
#'   \item{morekids}{indicator for having more than two children vs. exactly two children.}
#'   \item{samesex}{indicator for the first two children having the same sex (male-male or female-female)}
#'   \item{yob}{the year the woman was born}
#' }
#' @source Derived from Angrist and Evans (1998, The American Economic Review).
"AE"

#' ivmte Simulated Data
#'
#' @format A data frame with 5,000 rows and 14 columns.
#' \describe{
#'   \item{y}{binary outcome variable}
#'   \item{d}{binary treatment variable}
#'   \item{z}{instrument that takes the value 0, 1, 2, or 3}
#'   \item{x}{covariate x that takes integer values from 1 to 10}
#' }
#' @source Simulated --- see code in data/ivmteSimData.R.
"ivmteSimData"
