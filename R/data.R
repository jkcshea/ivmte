#' Mosquito data set
#'
#' A simulated population-level data set characterizing the effect of
#' purchasing mosquito nets, and the likelihood of contracting
#' malaria. A randomly assigned subsidy to purchase mosquito nets was
#' provided as part of an experiment.
#'
#' @format A data frame with 400 rows and 7 columns.
#' \describe{
#'   \item{i}{index for observation}
#'   \item{z}{categorical variable for level of subsidy}
#'   \item{pz}{probability of purchasing mosquito net}
#'   \item{d}{indicator for whether or not mosquito net was purchased}
#' 
#'   \item{ey0}{counterfactual probability of contracting malaria
#'   conditional on not purchasing a mosquito net.}
#'
#'   \item{ey1}{counterfactual probability of contracting malaria
#'   conditional on purchasing a mosquito net.}
#' 
#'   \item{ey}{the observed probability of contracting malaria.}
#' }
#' @source Simulated, based on Mogstad, Torgovitsky (2017).
"dtm"

#' Covariates data set
#'
#' A simulated population-level data set characterizing the effect of
#' a treatment on an outcome. The data includes a treatment indicator,
#' two covariates and two instruments.
#'
#' @format A data frame with 10,000 rows and 9 columns.
#' \describe{
#'   \item{x1}{covariate 1}
#'   \item{x2}{covariate 2}
#'   \item{z1}{instrument 1}
#'   \item{z2}{instrument 2}
#'   \item{p}{probability of purchasing mosquito net}
#'   \item{d}{indicator for treatment (d = 1) versus control (d = 0) group}
#'   \item{ey0}{counterfactual outcome when not a recipient of treatment}
#'   \item{ey1}{counterfactual outcome when a recipient of treatment}
#'   \item{ey}{the observed outcome}
#' }
#' @source Simulated.
"dtcf"
