gendist1 <- function(subN = 5, p1 = 0.4, p2 = 0.6) {
    ## subN <- 5
    ## p1 <- 0.4
    ## p2 <- 0.6

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
    plist <- polyparse.mst(~ u, data = dt, uname = u)
    plist0 <- gengamma.mst(plist,
                           lb = dt$p,
                           ub = 1,
                           multiplier = 1 / (1 - dt$p),
                           means = FALSE)
    plist1 <- gengamma.mst(plist,
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

gendist2 <- function(subN = 5, p1 = 0.4, p2 = 0.6, p3 = 0.8) {
    ## subN <- 5
    ## p1 <- 0.4
    ## p2 <- 0.6
    ## p3 <- 0.8

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
    plist <- polyparse.mst(~ u, data = dt, uname = u)
    plist0 <- gengamma.mst(plist,
                           lb = dt$p,
                           ub = 1,
                           multiplier = 1 / (1 - dt$p),
                           means = FALSE)
    plist1 <- gengamma.mst(plist,
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

gendist3 <- function(subN = 5, p1 = 0.4, p2 = 0.6) {

    dt <- data.frame(i  = rep(seq(1, subN), 2),
                     z  = rep(c(1, 2), each = subN),
                     p = rep(c(p1, p2), each = subN))

    ## Assign treatment
    dt$d <- 0
    dt[dt$z == 1 & dt$i <= (p1 * subN), "d"] <- 1
    dt[dt$z == 2 & dt$i <= (p2 * subN), "d"] <- 1

    ## Now construct Y as E[Y | X, D, Z]
    plist <- polyparse.mst(~ 0 + 1, data = dt, uname = u)
    plist0 <- gengamma.mst(plist,
                           lb = dt$p,
                           ub = 1,
                           multiplier = 1 / (1 - dt$p),
                           means = FALSE)
    plist1 <- gengamma.mst(plist,
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

