#' Estimating propensity scores
#'
#' This function estimates the propensity of taking up treatment. The
#' user can choose from fitting a linear probability model, a logit
#' model, or a probit model. The function can also be used to generate
#' a table of propensity scores for a given set of covariates and
#' excluded variables. This was incorporated to account for the LATE
#' being a target parameter. Specifically, if the argument
#' \code{formula} is the name of a variable in \code{data}, but the
#' target parameter is not the LATE, then no propensity model is
#' returned. If the target parameter is the LATE, then then the
#' propensity model is simply the empirical distribution of propensity
#' scores in the data conditioned on the set of covariates declared in
#' \code{late.X} and \code{late.Z}.
#'
#' @param formula Formula characterizing probability model. If a
#'     variable in the data already contains the propensity scores,
#'     input the variable as a one-sided formula. For example, if the
#'     variable \code{pz} contains the propensity score, input
#'     \code{formula = ~ pz}.
#' @param data \code{data.frame} with which to estimate the model.
#' @param link Link function with which to estimate probability
#'     model. Can be chosen from "linear", "logit", or "probit".
#' @param env environment, the environment for the original
#'     propensity score formula.
#' @return A vector of propensity scores for each observation, as well
#'     as a 'model'. If the user inputs a formula characterizing the
#'     model for taking up treatment, then the \code{lm}/\code{glm}
#'     object is returned. If the user declares a variable in the data
#'     set to be used as the propensity score, then a
#'     \code{data.frame} containing the propensity score for each
#'     value of the covariates in the probability model is returned.
#'
#' @examples
#' dtm <- ivmte:::gendistMosquito()
#'
#' ## Declaring a probability model.
#' propensity(formula = d ~ z,
#'                data = dtm,
#'                link = "linear")
#'
#' ## Declaring a variable to be used instead
#' propensity(formula = ~ pz,
#'                data = dtm,
#'                link = "linear")
#'
#' @export
propensity <- function(formula, data, link = "logit", env = parent.frame()) {
    formula <- Formula::as.Formula(formula)
    environment(formula) <- env
    propVars <- all.vars(formula)
    ## If two-sided formula is provided, estimate propensity score
    ## accordingly
    if (length(formula)[1] == 1 & length(formula)[2] == 1) {
        ## obtain design matrix
        if (link == "linear") prop <-  lm(formula, data, x = FALSE, y = FALSE)
        if (link == "logit")  prop <- glm(formula,
                                          family = binomial(link = "logit"),
                                          data, x = FALSE, y = FALSE)
        if (link == "probit") prop <- glm(formula,
                                          family = binomial(link = "probit"),
                                          data, x = FALSE, y = FALSE)
        phat <- prop$fitted.values
        phat[which(phat < 0)] <- 0
        phat[which(phat > 1)] <- 1
        return(list(model = prop, phat = phat))
    } else if (length(formula)[1] == 0 & length(formula)[2] == 1) {
        ## If one-sided formula is provided, containing only one
        ## variable, then that variable is the propensity score.
        ## Check that only one variable is submitted. This variable is
        ## the one in the data containing the propensity scores
        if (length(all.vars(formula)) > 1) {
            stop(gsub("\\s+", " ",
                      paste0("'formula' argument characterizing the
                      propensity score model must either be a two sided
                      formula (if the propensity score is to be estimated from
                      the data), or a one-sided formula containing a single
                      variable on the RHS (where the variable listed is included
                      in the data, and corresponds to propensity scores.")))
        }
        pname <- all.vars(formula)        
        ## Check that the variable is indeed a propensity,
        ## i.e. bounded between 0 and 1
        phat <- data[[pname]]
        names(phat) <- rownames(data)
        if (max(phat) > 1 | min(phat) < 0) {
            stop('Provided propensity scores are not between 0 and 1.')
        }
        model <- NULL
        return(list(variable = pname, model = model, phat = phat))
    } else {
        stop(gsub("\\s+", " ",
                  paste0("'formula' argument characterizing the
                  propensity score model must either be a two sided formula
                  (if the propensity score is to be estimated from the data),
                  or a one-sided formula containing a single variable on the RHS
                  (where the variable listed is included in the data, and
                  corresponds to propensity scores.")))
    }
}
