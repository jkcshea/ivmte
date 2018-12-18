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
#' @param pi1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param pi2 the probability of treatment for those with the
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
    plist0 <- gengamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- gengamma(plist,
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
#' @param pi1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param pi2 the probability of treatment for those with the
#'     instrument Z = 2.
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
#' @param pi1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param pi2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @param p13 the probability of treatment for those with the
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
    plist0 <- gengamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- gengamma(plist,
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
#' @param pi1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param pi2 the probability of treatment for those with the
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
    plist0 <- gengamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- gengamma(plist,
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
#' @param pi1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param pi2 the probability of treatment for those with the
#'     instrument Z = 2.
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
#' @param pi1 the probability of treatment for those with the
#'     instrument Z = 1.
#' @param pi2 the probability of treatment for those with the
#'     instrument Z = 2.
#' @param p13 the probability of treatment for those with the
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
    plist0 <- gengamma(plist,
                       lb = dt$p,
                       ub = 1,
                       multiplier = 1 / (1 - dt$p),
                       means = FALSE)
    plist1 <- gengamma(plist,
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


