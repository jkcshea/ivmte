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
#' @param formula Formula characterizing probability model, or name of
#'     variable in data set containing propensity scores.
#' @param data \code{data.frame} with which to estimate the model.
#' @param link Link function with which to estimate probability
#'     model. Can be chosen from "linear", "logit", or "probit".
#' @param late.Z A vector of variable names of excluded
#'     variables. This is required when the target parameter is the
#'     LATE.
#' @param late.X A vector of variable names of non-excluded
#'     variables. This is required when the target parameter is the
#'     LATE, and the estimation procedure will condition on these
#'     variables.
#' @return A vector of propensity scores for each observation, as well
#'     as a `model'. If the user inputs a formula characterizing the
#'     model for taking up treatment, then the \code{lm}/\code{glm}
#'     object is returned. If the user declares a variable in the data
#'     set to be used as the propensity score, then a
#'     \code{data.frame} containing the propensity score for each
#'     value of the covariates in the probability model is returned.
#'
#' @examples
#' Declaring a probability model.
#' propensity.mst(formula = d ~ x1 + x2 + z1 + z2,
#'                data = data,
#'                link = "logit")
#'
#' Declaring a variable to be used instead
#' propensity.mst(formula = p,
#'                data = data,
#'                late.Z = c(z1, z2),
#'                late.X = c(x1, x2))
#'
#' @export
propensity.mst <- function(formula, data, link = "linear", late.Z,
                           late.X) {

    fname <- deparse(substitute(formula))

    ## If formula is provided, estimate propensity score accordingly
    if (class_formula(as.formula(fname))) {
        ## obtain design matrix
        if (link == "linear") prop <-  lm(formula, data)
        if (link == "logit")  prop <- glm(formula,
                                          family = binomial(link = "logit"),
                                          data)
        if (link == "probit") prop <- glm(formula,
                                          family = binomial(link = "probit"),
                                          data)

        phat <- prop$fitted.values

        phat[which(phat < 0)] <- 0
        phat[which(phat > 1)] <- 1
        return(list(model = prop, phat = phat))
    } else {
        ## If variable name is provided, check if it is indeed a
        ## propensity score
        phat <- data[[fname]]
        names(phat) <- rownames(data)

        if (max(phat) > 1 | min(phat) < 0) {
            stop('Provided propensity scores are not between 0 and 1.')
        }
        model <- NULL
        if (hasArg(late.Z)) {
            ## Generate a matrix of propensity scores corresponding to
            ## each X and Z---this only works if the X and Z are
            ## declared for the LATE variables.
            late.Z <- restring(substitute(late.Z), substitute = FALSE)

            if(hasArg(late.X))  late.X <- restring(substitute(late.X),
                                                   substitute = FALSE)
            if(!hasArg(late.X)) late.X <- NULL
            model <- unique(data[, c(late.X, late.Z, fname)])

            ## There should be no duplicate (X, Z)s in the model
            ## (i.e. for each set of (X, Z), there should only be one
            ## propensity score)
            if(max(duplicated(model[, c(late.X, late.Z)]))) {
                stop(gsub("\\s+", " ",
                          "Custom propensity score is not well defined: some
                          (X, Z) values are mapped to more than one
                          propensity."))
           }
        }
        return(list(model = model, phat = phat))
    }
}
