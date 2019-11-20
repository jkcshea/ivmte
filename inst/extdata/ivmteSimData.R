set.seed(1)
n <- 5000
u <- runif(n)
z <- rbinom(n, 3, .5)
x <- as.numeric(cut(rnorm(n), 10)) # normal discretized into 10 bins
d <- as.numeric(u < z*.25 + .01*x)

v0 <- rnorm(n) + .2*u
m0 <- 0
y0 <- as.numeric(m0 + v0 + .1*x > 0)

v1 <- rnorm(n) - .2*u
m1 <- .5
y1 <- as.numeric(m1 + v1 - .3*x > 0)

y <- d*y1 + (1-d)*y0

ivmteSimData <- data.frame(y,d,z,x)

library("devtools")
use_data(ivmteSimData, overwrite = TRUE)