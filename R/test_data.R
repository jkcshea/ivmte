#' Generate test distribution 1
#'
#' This function generates a data set for testing purposes. There is a
#' single instrument that takes on values of 1 or 2, and the
#' distribution of the values for the binary instrument is
#' uniform. The MTRs are m0 ~ 0 + u and m1 ~ 1 + u. All unobservables
#' u are integrated out.
#' @param subN integer, default set to 5. This is the number of
#'     individuals possessing each value of the instrument. So the
#'     total number of observations is subN * 2.
#' @param p1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param p2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @return a data.frame.
gendist1 <- function(subN = 5, p1 = 0.4, p2 = 0.6) {
    dt <- data.frame(i  = rep(seq(1, subN), 2),
                     z  = rep(c(1, 2), each = subN),
                     p = rep(c(p1, p2), each = subN))

    ## Assign treatment
    dt$d <- 0
    dt[dt$z == 1 & dt$i <= (p1 * subN), "d"] <- 1
    dt[dt$z == 2 & dt$i <= (p2 * subN), "d"] <- 1

    ## Now construct Y as E[Y | X, D, Z]
    ## m0 = 1 + u
    ## m1 = 1 + u
    plist <- polyparse(~ u, data = dt, uname = u)
    plist0 <- genGamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- genGamma(plist,
                       lb = 0,
                       ub = dt$p,
                       multiplier = 1 / dt$p,
                       means = FALSE)

    m0coef <- c(0, 6)
    m1coef <- c(7, 8)

    dt$ey0 <- plist0 %*% m0coef
    dt$ey1 <- plist1 %*% m1coef
    dt$ey  <-  dt$d * dt$ey1 + (1 - dt$d) * dt$ey0

    ## Convert data back into data.frame
    return(dt)
}

#' Generate test distribution 1 with errors
#'
#' This function generates a data set for testing purposes. There is a
#' single instrument that takes on values of 1 or 2, and the
#' distribution of the values for the binary instrument is
#' uniform. The MTRs are m0 ~ 0 + u and m1 ~ 1 + u.
#'@param N integer, default set to 100. Total number of observations
#'     in the data.
#' @param subN , default set to 0.5. This is the probability the agent
#'     will have Z = 1.
#' @param p1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param p2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @param v0.sd numeric, standard deviation of error term for
#'     counterfactual \code{D = 0}
#' @param v1.sd numeric, standard deviation of error term for
#'     counterfactual \code{D = 1}
#' @return a data.frame.
gendist1e <- function(N = 100, subN = 0.5, p1 = 0.4, p2 = 0.6, v0.sd = 0.5,
                      v1.sd = 0.75) {

    dt <- data.frame(z = runif(N),
                     u = runif(N))

    ## First construct the instruments independently from U
    dt$z  <- as.integer(dt$z > subN) + 1

    dt$d <- 0
    dt[dt$z == 1 & dt$u <= p1, "d"] <- 1
    dt[dt$z == 2 & dt$u <= p2, "d"] <- 1

    ## Now construct the error terms. The error terms are normally
    ## distributed, centered at the draw of u if D = 0, and -u if D =
    ## 1.
    dt$v0 <- sapply(dt[, "u"], rnorm, n = 1, sd = v0.sd)
    dt$v1 <- sapply(-dt[, "u"], rnorm, n = 1, sd = v1.sd)

    ## Now construct the Y values
    ## m0 = 0 + u
    ## m1 = 1 + u
    dt[dt$d == 0, "y"] <- 0 + 6 * dt[dt$d == 0, "u"] + dt[dt$d == 0, "v0"]
    dt[dt$d == 1, "y"] <- 7 + 8 * dt[dt$d == 1, "u"] + dt[dt$d == 1, "v1"]

    ## Label the types
    dt[dt$d == 1 & dt$z == 1, "type"] <- 1
    dt[dt$d == 0 & dt$z == 1, "type"] <- 2
    dt[dt$d == 1 & dt$z == 2, "type"] <- 3
    dt[dt$d == 0 & dt$z == 2, "type"] <- 4

    return(dt)
}

#' Generate test distribution 2
#'
#' This function generates a data set for testing purposes. There is a
#' single instrument that takes on values of 1, 2, or 3, and the
#' distribution of the values for the binary instrument is
#' uniform. The MTRs are m0 ~ 1 + u and m1 ~ 1 + u. All unobservables
#' u are integrated out.
#' @param subN integer, default set to 5. This is the number of
#'     individuals possessing each value of the instrument. So the
#'     total number of observations is subN * 2.
#' @param p1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param p2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @param p3 the probability of treatment for those with the
#'     instrument Z = 3.
#' @return a data.frame.
gendist2 <- function(subN = 5, p1 = 0.4, p2 = 0.6, p3 = 0.8) {

    dt <- data.frame(i  = rep(seq(1, subN), 3),
                     z  = rep(c(1, 2, 3), each = subN),
                     p = rep(c(p1, p2, p3), each = subN))

    ## Assign treatment
    dt$d <- 0
    dt[dt$z == 1 & dt$i <= p1 * subN, "d"] <- 1
    dt[dt$z == 2 & dt$i <= p2 * subN, "d"] <- 1
    dt[dt$z == 3 & dt$i <= p3 * subN, "d"] <- 1

    ## Now construct Y as E[Y | X, D, Z]
    ## m0 = 1 + u
    ## m1 = 1 + u
    plist <- polyparse(~ u, data = dt, uname = u)
    plist0 <- genGamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- genGamma(plist,
                       lb = 0,
                       ub = dt$p,
                       multiplier = 1 / dt$p,
                       means = FALSE)

    m0coef <- c(0, 6)
    m1coef <- c(7, 8)

    dt$ey0 <- plist0 %*% m0coef
    dt$ey1 <- plist1 %*% m1coef
    dt$ey  <- dt$d * dt$ey1 + (1 - dt$d) * dt$ey0

    ## Convert data back into data.frame
    return(dt)
}


#' Generate test distribution 3
#'
#' This function generates a data set for testing purposes. There is a
#' single instrument that takes on values of 1 and 2, and the
#' distribution of the values for the binary instrument is
#' uniform. The MTRs are m0 ~ 1 and m1 ~ 1. All unobservables u are
#' integrated out.
#' @param subN integer, default set to 5. This is the number of
#'     individuals possessing each value of the instrument. So the
#'     total number of observations is subN * 2.
#' @param p1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param p2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @return a data.frame.
gendist3 <- function(subN = 5, p1 = 0.4, p2 = 0.6) {

    dt <- data.frame(i  = rep(seq(1, subN), 2),
                     z  = rep(c(1, 2), each = subN),
                     p = rep(c(p1, p2), each = subN))

    ## Assign treatment
    dt$d <- 0
    dt[dt$z == 1 & dt$i <= (p1 * subN), "d"] <- 1
    dt[dt$z == 2 & dt$i <= (p2 * subN), "d"] <- 1

    ## Now construct Y as E[Y | X, D, Z]
    plist <- polyparse(~ 0 + 1, data = dt, uname = u)
    plist0 <- genGamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- genGamma(plist,
                       lb = 0,
                       ub = dt$p,
                       multiplier = 1 / dt$p,
                       means = FALSE)

    m0coef <- c(3)
    m1coef <- c(9)

    dt$ey0 <- plist0 %*% m0coef
    dt$ey1 <- plist1 %*% m1coef
    dt$ey  <-  dt$d * dt$ey1 + (1 - dt$d) * dt$ey0

    ## Convert data back into data.frame
    return(dt)
}


#' Generate test distribution 3 with errors
#'
#' This function generates a data set for testing purposes. There is a
#' single instrument that takes on values of 1 or 2, and the
#' distribution of the values for the binary instrument is
#' uniform. The MTRs are m0 ~ 0 + u and m1 ~ 1 + u.
#'@param N integer, default set to 100. Total number of observations
#'     in the data.
#' @param subN , default set to 0.5. This is the probability the agent
#'     will have Z = 1.
#' @param p1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param p2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @param v0.sd numeric, standard deviation of error term for
#'     counterfactual \code{D = 0}
#' @param v1.sd numeric, standard deviation of error term for
#'     counterfactual \code{D = 1}
#' @return a data.frame.
gendist3e <- function(N = 100, subN = 0.5, p1 = 0.4, p2 = 0.6, v0.sd = 0.5,
                      v1.sd = 0.75) {

    dt <- data.frame(z = runif(N),
                     u = runif(N))

    ## First construct the instruments independently from U
    dt$z  <- as.integer(dt$z > subN) + 1

    dt$d <- 0
    dt[dt$z == 1 & dt$u <= p1, "d"] <- 1
    dt[dt$z == 2 & dt$u <= p2, "d"] <- 1

    ## Now construct the error terms. The error terms are normally
    ## distributed, centered at the draw of u if D = 0, and -u if D =
    ## 1.
    dt$v0 <- sapply(dt[, "u"], rnorm, n = 1, sd = v0.sd)
    dt$v1 <- sapply(-dt[, "u"], rnorm, n = 1, sd = v1.sd)

    ## Now construct the Y values
    ## m0 = 3
    ## m1 = 9
    dt[dt$d == 0, "y"] <- 3 + dt[dt$d == 0, "v0"]
    dt[dt$d == 1, "y"] <- 9 + dt[dt$d == 1, "v1"]

    ## Label the types
    dt[dt$d == 1 & dt$z == 1, "type"] <- 1
    dt[dt$d == 0 & dt$z == 1, "type"] <- 2
    dt[dt$d == 1 & dt$z == 2, "type"] <- 3
    dt[dt$d == 0 & dt$z == 2, "type"] <- 4

    return(dt)
}


#' Generate test distribution 4
#'
#' This function generates a data set for testing purposes. There is a
#' single instrument that takes on values of 1, 2, and 3, and the
#' distribution of the values for the binary instrument is
#' uniform. The MTRs are m0 ~ 1 and m1 ~ 1. All unobservables u are
#' integrated out.
#' @param subN integer, default set to 5. This is the number of
#'     individuals possessing each value of the instrument. So the
#'     total number of observations is subN * 2.
#' @param p1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param p2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @param p3 the probability of treatment for those with the
#'     instrument Z = 3.
#' @return a data.frame.
gendist4 <- function(subN = 5, p1 = 0.4, p2 = 0.6, p3 = 0.8) {

    dt <- data.frame(i  = rep(seq(1, subN), 3),
                     z  = rep(c(1, 2, 3), each = subN),
                     p = rep(c(p1, p2, p3), each = subN))

    ## Assign treatment
    dt$d <- 0
    dt[dt$z == 1 & dt$i <= (p1 * subN), "d"] <- 1
    dt[dt$z == 2 & dt$i <= (p2 * subN), "d"] <- 1
    dt[dt$z == 3 & dt$i <= (p3 * subN), "d"] <- 1

    ## Now construct Y as E[Y | X, D, Z]
    plist <- polyparse(~ 0 + 1, data = dt, uname = u)
    plist0 <- genGamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- genGamma(plist,
                       lb = 0,
                       ub = dt$p,
                       multiplier = 1 / dt$p,
                       means = FALSE)

    m0coef <- c(3)
    m1coef <- c(9)

    dt$ey0 <- plist0 %*% m0coef
    dt$ey1 <- plist1 %*% m1coef
    dt$ey  <-  dt$d * dt$ey1 + (1 - dt$d) * dt$ey0

    ## Convert data back into data.frame
    return(dt)
}



#' Generate mosquito data set
#'
#' This code generates the population level data in Mogstad, Santos,
#' Torgovitsky (2018), i.e. the mosquito data set used as the running
#' example.
#' 
#' @return data.frame.
gendist_mosquito <- function() {

    set.seed(1)

    subN <- 100

    p1 <- 0.12
    p2 <- 0.29
    p3 <- 0.48
    p4 <- 0.78

    dt <- data.frame(i  = rep(seq(1, subN), 4),
                     z  = rep(c(1, 2, 3, 4), each = subN),
                     pz = rep(c(p1, p2, p3, p4), each = subN))
    
    ## Assign treatment
    dt$d <- 0
    dt[dt$z == 1 & dt$i <= 12, "d"] <- 1
    dt[dt$z == 2 & dt$i <= 29, "d"] <- 1
    dt[dt$z == 3 & dt$i <= 48, "d"] <- 1
    dt[dt$z == 4 & dt$i <= 78, "d"] <- 1

    ## Now construct Y as E[Y | X, D, Z]
    plist <- polyparse(~ u + I(u ^ 2), data = dt, uname = u)
    plist0 <- genGamma(plist,
                           lb = dt$pz,
                           ub = 1,
                           multiplier = 1 / (1 - dt$pz),
                           means = FALSE)
    plist1 <- genGamma(plist,
                           lb = 0,
                           ub = dt$pz,
                           multiplier = 1 / dt$pz,
                           means = FALSE)

    m0coef <- c(0.9, -1.1, 0.3)
    m1coef <- c(0.35, -0.3, -0.05)

    dt$ey0 <- plist0 %*% m0coef
    dt$ey1 <- plist1 %*% m1coef
    dt$ey  <- dt$d * dt$ey1 + (1 - dt$d) * dt$ey0

    ## Convert data back into data.frame
    dtm <- data.frame(dt)

    return(dtm)
}

#' Generate test data set with covariates
#' 
#' This code generates population level data to test the estimation
#' function. This data includes covariates. The data generated will
#' have already integrated over the unobservable terms U, where U | X,
#' Z ~ Unif[0, 1].
#' 
#' @return a list of two data.frame objects. One is the distribution
#'     of the simulated data, the other is the full simulated data
#'     set.
gendist_covariates <- function() {

    set.seed(1)

    supp_x1 <- c(0, 1)
    supp_x2 <- c(1, 2, 3)
    supp_z1 <- c(0, 1)
    supp_z2 <- c(1, 2, 3)
    
    dtc <- data.frame(expand.grid(supp_x1, supp_x2,
                                  supp_z1, supp_z2))
    
    colnames(dtc) <- c("x1", "x2", "z1", "z2")

    ## Generate propensity scores
    g0 <- -0.5
    g1 <- -0.4
    g2 <- 0.2
    g3 <- -1
    g4 <- 0.3
    
    dtc$latent <- g0 + g1 * dtc$x1 + g2 * dtc$x2 +
        g3 * dtc$z1 + g4 * dtc$z2
    dtc$p <- round(1 / (1 + exp(-dtc$latent)), 2)
    
    ## Generate the counterfactual outcomes Note: the code is prepared
    ## such that the polynomials are all in terms of u. Thus, terms which
    ## have u components with the same power are all grouped together.
    ##
    ## So I use the following specifications
    ## m0 = 0.3 + 0.4 * x1 - 0.1 * x2 * u - 0.2 * x2 * u^2
    ## m1 = 0.5 + 0.2 * x1 - 0.1 * x1 * x2 - 0.02 * u +
    ##        0.3 * x1 * u - 0.05 * x2 * u^2

    plist0 <- polyparse(~ x1 + I(x2 * u) + I(x2 * u^2),
                            data = dtc,
                            uname = u)
    plist1 <- polyparse(~ x1 + I(x1 * x2) + u + I(x1 * u) + I(x2 * u^2),
                            data = dtc,
                            uname = u)

    glist0 <- genGamma(plist0,
                           lb = dtc$p,
                           ub = 1,
                           multiplier = 1 / (1 - dtc$p),
                           means = FALSE)
    glist1 <- genGamma(plist1,
                           lb = 0,
                           ub = dtc$p,
                           multiplier = 1 / dtc$p,
                           means = FALSE)

    g0coef <- c(0.3, 0.4, -0.1, -0.2)
    g1coef <- c(0.5, 0.2, -0.1, -0.02, 0.3, -0.05)
   
    dtc$ey0 <- glist0 %*% g0coef
    dtc$ey1 <- glist1 %*% g1coef
    
    ## Generate distribution

    ## I want to generate the data to allow for correlation across the
    ## variables. But then designing the distribution manually is too
    ## difficult. So instead, I try an approach where I build the
    ## distributions marginally, allowing for correlations to occur, and
    ## then normalize the distribution.

    dtc$f <- 0
  
    dtc[dtc$x1 == 0, "f"] <- dtc[dtc$x1 == 0, "f"] + 0.1
    dtc[dtc$x1 == 1, "f"] <- dtc[dtc$x1 == 1, "f"] + 0.13
    
    dtc[dtc$x2 == 1, "f"] <- dtc[dtc$x2 == 1, "f"] + 0.05
    dtc[dtc$x2 == 2, "f"] <- dtc[dtc$x2 == 2, "f"] + 0.1
    dtc[dtc$x2 == 3, "f"] <- dtc[dtc$x2 == 3, "f"] + 0.01
    
    dtc[dtc$x1 == 0 & dtc$z1 == 1, "f"] <-
        dtc[dtc$x1 == 0 & dtc$z1 == 1, "f"] - 0.03 
    dtc[dtc$x1 == 0 & dtc$z2 == 2, "f"] <-
        dtc[dtc$x1 == 0 & dtc$z2 == 2, "f"] - 0.01
    dtc[dtc$x1 == 0 & dtc$z2 == 3, "f"] <-
        dtc[dtc$x1 == 0 & dtc$z2 == 3, "f"] - 0.02

    dtc[dtc$z1 == 1 & dtc$z2 == 2, "f"] <-
        dtc[dtc$z1 == 1 & dtc$z2 == 2, "f"] + 0.02
    dtc[dtc$z1 == 1 & dtc$z2 == 3, "f"] <-
        dtc[dtc$z1 == 1 & dtc$z2 == 3, "f"] + 0.01

    dtc[dtc$z2 == 2, "f"] <- dtc[dtc$z2 == 2, "f"] + 0.05
    dtc[dtc$z2 == 3, "f"] <- dtc[dtc$z2 == 3, "f"] + 0.01
    dtc[dtc$x2 == 2 & dtc$z2 == 2, "f"] <-
        dtc[dtc$x2 == 2 & dtc$z2 == 2, "f"] + 0.01
    dtc[dtc$x2 == 2 & dtc$z2 == 3, "f"] <-
        dtc[dtc$x2 == 2 & dtc$z2 == 3, "f"] + 0.01
    dtc[dtc$x2 == 3 & dtc$z2 == 2, "f"] <-
        dtc[dtc$x2 == 3 & dtc$z2 == 2, "f"] + 0.02
    dtc[dtc$x2 == 3 & dtc$z2 == 3, "f"] <-
        dtc[dtc$x2 == 3 & dtc$z2 == 3, "f"] + 0.04

    dtc$f <- dtc$f / sum(dtc$f)
    
    ## Check distribution (requires data.table)
    ## length(unique(dtc$f))
    ## dtc[, sum(f), by = x1]
    ## dtc[, sum(f), by = x2]
    ## dtc[, sum(f), by = z1]
    ## dtc[, sum(f), by = z2]
    
    ## Since distributions were arbitrarily designed anyways, I will round
    ## them to make it easier to work with. I will make sure that they sum
    ## to 1.
    dtc$f <- round(dtc$f, 2)
    sum(dtc$f) ## Sum of rounded densities is 1.02. Adjust the
               ## densities so they sum to 1.

    dtc[dtc$x1 == 1 & dtc$x2 == 1 & dtc$z1 == 0 & dtc$z2 == 1, "f"] <-
        dtc[dtc$x1 == 1 & dtc$x2 == 1 & dtc$z1 == 0 & dtc$z2 == 1, "f"] - 0.01
    dtc[dtc$x1 == 1 & dtc$x2 == 3 & dtc$z1 == 1 & dtc$z2 == 3, "f"] <-
        dtc[dtc$x1 == 1 & dtc$x2 == 3 & dtc$z1 == 1 & dtc$z2 == 3, "f"] - 0.01
    dtc$f <- round(dtc$f, 2)
    
    ## Expand data set according to distribution

    ## We multiply each row by 100 * (100 * f). The first 100 is so that
    ## we assign the probability of treatment and control according to the
    ## variable p. The (100 * f) is so that we can get the distribution of
    ## covariates to be according to f.

    ## dtc[, multiplier := 100 * 100 * f]
    ## sum(dtc$multiplier)
    ## multipler <- dtc$multiplier
    ## dtcf <- dtc[rep(seq_len(nrow(dtc)), multiplier),]

    dtc[, "multiplier"] <- 100 * 100 * dtc$f
    sum(dtc$multiplier)
    dtcf <- dtc[rep(seq_len(nrow(dtc)), dtc$multiplier), ]

    ## Check if the distribution of this data set matches that of what we
    ## generated above (requires data.table)
    ## round(dtcf[, .N/nrow(dtcf), by = x1], 2) == round(dtc[, sum(f), by = x1], 2)
    ## round(dtcf[, .N/nrow(dtcf), by = x2], 2) == round(dtc[, sum(f), by = x2], 2)
    ## round(dtcf[, .N/nrow(dtcf), by = z1], 2) == round(dtc[, sum(f), by = z1], 2)
    ## round(dtcf[, .N/nrow(dtcf), by = z2], 2) == round(dtc[, sum(f), by = z2], 2)

    ## Now assign treatment
    dtcf$d <- 0
    for (ix1 in supp_x1) {
        for (ix2 in supp_x2) {
            for (iz1 in supp_z1) {
                for (iz2 in supp_z2) {
                    N <- nrow(dtcf[dtcf$x1 == ix1 &
                                   dtcf$x2 == ix2 &
                                   dtcf$z1 == iz1 &
                                   dtcf$z2 == iz2, ])
                    dtcf[dtcf$x1 == ix1 &
                         dtcf$x2 == ix2 &
                         dtcf$z1 == iz1 &
                         dtcf$z2 == iz2, "i"] <- seq(1, N)
                    dtcf[dtcf$x1 == ix1 &
                         dtcf$x2 == ix2 &
                         dtcf$z1 == iz1 &
                         dtcf$z2 == iz2, "dcut"] <-
                        round(dtcf[dtcf$x1 == ix1 &
                                   dtcf$x2 == ix2 &
                                   dtcf$z1 == iz1 &
                                   dtcf$z2 == iz2, "p"] * N)
                }
            }
        }
    }
    dtcf[dtcf$i <= dtcf$dcut, "d"] <- 1
    
    ## Check if empirical distribution differs from population
    ## (requires data.table)

    ## failed <- which(dtcf[, mean(d), by = .(x1, x2, z1, z2)]$V1 !=
    ##                 dtcf[, mean(p), by = .(x1, x2, z1, z2)]$V1)
    ## print(failed)

    ## We want this list to be empty (i.e. a failure occurs when the
    ## empirical distribution of treatment within group differs from
    ## what is specified in the population).

    ## Assign outcomes
    dtcf$ey <- dtcf$d * dtcf$ey1 + (1 - dtcf$d) * dtcf$ey0

    return(list(data.full = dtcf,
                data.dist = dtc))
}

#' Generate test data set with splines
#' 
#' This code generates population level data to test the
#' estimation function. This data set incorporates splines in the
#' MTRs.
#' 
#' The distribution of the data is as follows
#'
#'        |     Z
#'  X/Z   |  0     1
#' _______|___________
#'     -1 | 0.1   0.1
#'        |
#'  X   0 | 0.2   0.2
#'        |
#'      1 | 0.1   0.2
#'
#' The data presented below will have already integrated over the
#' unobservable terms U, and U | X, Z ~ Unif[0, 1].
#'
#' The propensity scores are generated according to the model
#'
#' p(x, z) = 0.5 - 0.1 * x + 0.2 * z
#'
#'        |     Z
#' p(X,Z) |  0     1
#' _______|___________
#'     -1 | 0.6   0.8
#'        |
#'  X   0 | 0.5   0.7
#'        |
#'      1 | 0.4   0.6
#'
#' The lowest common multiple of the first table is 12. The lowest
#' common multiple of the second table is 84. It turns out that 840 *
#' 5 = 4200 observations is enough to generate the population data
#' set, such that each group has a whole-number of observations.
#'
#' The MTRs are defined as follows:
#'
#' y1 ~ beta0 + beta1 * x + uSpline(degree = 2,
#'                                 knots  = c(0.3, 0.6),
#'                                 intercept = FALSE)
#'
#' The coefficients (beta1, beta2), and the coefficients on the
#' splines, will be defined below.
#'
#' y0 = x : uSpline(degree = 0,
#'                  knots  = c(0.2, 0.5, 0.8),
#'                  intercept = TRUE)
#'      + uSpline(degree = 1,
#'                knots  = c(0.4),
#'                intercept = TRUE)
#'      + beta3 * I(u ^ 2)
#'
#' The coefficient beta3, and the coefficients on the splines, will be
#' defined below.
#'
#' @return a list of two data.frame objects. One is the distribution
#'     of the simulated data, the other is the full simulated data
#'     set.
gendist_splines <- function() {

    set.seed(10L)

    ## Declare coefficients
    pCoef    <- c(0.5, -0.1, 0.2) ## Propensity score coefficients
    y1Coef   <- c(30, 35) ## Beta coefficients for y1
    y0Coef   <- c(20) ## Beta coefficients for y0
    u1s1Coef <- c(-75, 45, 75, 55) ## Coefficients for spline in y1
    u0s1Coef <- c(25, 60, 45, 30) ## Coefficients for first spline in y0
    u0s2Coef <- c(65, 70, 40) ## Coefficients for second spline in y0

    ## Generate distribution
    distr <- data.frame(group = seq(1, 6),
                        x = rep(c(-1, 0, 1), times = 2),
                        z = rep(c(0, 1), each = 3),
                        f = c(0.1, 0.2, 0.1, 0.1, 0.3, 0.2))

    distr$p <- pCoef[1] + pCoef[2] * distr$x + pCoef[3] * distr$z

    distr$ey1 <- y1Coef[1] * distr$p +
        y1Coef[2] * distr$x * distr$p +
        t(sapply(X = distr$p,
                 FUN = splineInt,
                 lb = 0,
                 degree = 2,
                 knots = c(0.3, 0.6),
                 intercept = FALSE)) %*% u1s1Coef

    distr$ey0 <- distr$x * t(sapply(X = distr$p,
                                    FUN = splineInt,
                                    ub = 1,
                                    degree = 0,
                                    knots = c(0.2, 0.5, 0.8),
                                    intercept = TRUE)) %*% u0s1Coef +
                           t(sapply(X = distr$p,
                                    FUN = splineInt,
                                    ub = 1,
                                    degree = 1,
                                    knots = c(0.4),
                                    intercept = TRUE)) %*% u0s2Coef +
                           y0Coef[1] / 3 * (1 - distr$p ^ 3)

    ## Expand distribution into data set
    N <- 4200
    distr$multiplier <- N * distr$f
    distr$controls <- distr$multiplier * (1 - distr$p)
    dtsf <- distr[rep(seq_len(nrow(distr)), distr$multiplier), c(1:7, 9)]
    rownames(dtsf) <- NULL

    ## Construct observed outcome
    dtsf$d <- 1
    dtsf$count <- 0
    for (g in unique(distr$group)) {
        dtsf[dtsf$group == g, "count"] <- seq(1, nrow(dtsf[dtsf$group == g, ]))
    }
    dtsf[round(dtsf$count) <= round(dtsf$controls), "d"] <- 0
    dtsf$ey <- (1 - dtsf$d) * dtsf$ey0 + dtsf$d * dtsf$ey1

    ## Check if the distribution is as planned
    for (z in c(0, 1)) {
        for (x in c(-1, 0, 1)) {
            a <- nrow(dtsf[dtsf$x == x & dtsf$z == z & dtsf$d == 1, ])
            b <- nrow(dtsf[dtsf$x == x & dtsf$z == z, ])
            cat("X:", x, ", Z:", z, ", a:", a, ", b:", b,
                ", Prob of treat:", a / b, "\n")
        }
    }

    ## Check if treatment effects are consistent
    as.numeric(t(distr$ey1 - distr$ey0) %*% distr$f)
    round(as.numeric(t(distr$ey1 - distr$ey0) %*% distr$f), 7) ==
        round(mean(dtsf$ey1 - dtsf$ey0), 7)
    mean(dtsf[dtsf$d == 1, "ey"]) - mean(dtsf[dtsf$d == 0, "ey"])

    ## Drop unnecssary columns
    dtsf[, c("group", "controls", "count")] <- NULL

    return(list(data.full = dtsf,
                data.dist = distr))
}

#' Generate basic data set for testing
#'
#' This code generates population level data to test the estimation
#' function. This is a simpler dataset, one in which we can more
#' easily estimate a correctly specified model.  The data presented
#' below will have already integrated over the # unobservable terms
#' U, where U | X, Z ~ Unif[0, 1].
#' 
#' @return a list of two data.frame objects. One is the distribution
#'     of the simulated data, the other is the full simulated data
#'     set.
gendist_basic <- function() {

    set.seed(10L)

    supp_x <- c(1, 2, 3)
    supp_z <- c(1, 2, 3, 4)

    dtb <- data.frame(expand.grid(supp_x, supp_z))
    colnames(dtb) <- c("x", "z")

    ## Generate propensity scores
    g0 <- 0.2
    g1 <- -0.10
    g2 <- 0.2
    dtb$p <- g0 + g1 * dtb$x + g2 * dtb$z

    ## Generate the counterfactual outcomes.
    ## m0 = 2 + x + 2 * u
    ## m1 = 6 + 5 * u^2

    plist0 <- polyparse(~ x + u,
                            data = dtb,
                            uname = u)
    plist1 <- polyparse(~ I(u^2),
                            data = dtb,
                            uname = u)

    ## Remember that you are generating E[m | D = 0] and E[m | D = 1],
    ## which are analogous to E[m | u > p(X, Z)] and E[m | u < p(X,
    ## Z)]. This is why you include those multipliers: You're
    ## integrating with respect to a conditional distribution.
    glist0 <- genGamma(plist0,
                           lb = dtb$p,
                           ub = 1,
                           multiplier = 1 / (1 - dtb$p),
                           means = FALSE)
    glist1 <- genGamma(plist1,
                           lb = 0,
                           ub = dtb$p,
                           multiplier = 1 / dtb$p,
                           means = FALSE)

    g0coef <- c(2, 1, 2)
    g1coef <- c(6, 5)

    dtb$ey0 <- glist0 %*% g0coef
    dtb$ey1 <- glist1 %*% g1coef   

    ## Generate distribution and expand data set
    p <- matrix(runif(12), ncol = 3)
    p <- p / sum(p)

    pmat <- matrix(c(0.02, 0.07, 0.02, 0.19,
                     0.09, 0.10, 0.17, 0.15,
                     0.01, 0.01, 0.16, 0.01),
                   ncol = 4)

    dtb$f <- c(pmat)
    dtb$multiplier <- dtb$f * 1000

    ## Check treatment effect (requires data.table)
    ## sum(dtb$ey1 * dtb$f) - sum(dtb$ey0 * dtb$f)
    ## sum(dtb[dtb$ey0 > dtb$ey1, f])

    ## Assign treatment
    dtbf <- dtb[rep(seq_len(nrow(dtb)), dtb$multiplier), ]
    dtbf$d <- 0
    for (ix in supp_x) {
        for (iz in supp_z) {
            N <- nrow(dtbf[dtbf$x == ix &
                           dtbf$z == iz, ])
            dtbf[dtbf$x == ix &
                 dtbf$z == iz, "i"] <- seq(1, N)
            dtbf[dtbf$x == ix &
                 dtbf$z == iz, "dcut"] <-
                round(dtbf[dtbf$x == ix &
                           dtbf$z == iz, "p"] * N)
        }
    }
    dtbf[dtbf$i <= dtbf$dcut, "d"] <- 1
    dtbf$dcut <- NULL

    ## Check if empirical distribution differs from population
    ## failed <- which(round(dtbf[, mean(d), by = .(x, z)]$V1, 2) !=
    ##                 round(dtbf[, mean(p), by = .(x, z)]$V1, 2))
    ## print(failed)

    ## Assign outcomes
    dtbf$ey <- dtbf$d * dtbf$ey1 + (1 - dtbf$d) * dtbf$ey0

    return(list(data.full = dtbf,
                data.dist = dtb))
}
