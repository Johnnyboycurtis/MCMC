
##---------------------------------------------------------------------------##
## Stat 591                                                                  ##
## Author: Ryan Martin (rgmartin@uic.edu)                                    ##
## R code solution to Gelman et al (2004), Exercise 3.12, page 99            ##
## Employs a very basic but clever strategy for approx posterior sampling    ##
##---------------------------------------------------------------------------##

library(MASS)


# Data from Gelman et al (2004), Table 2.2, page 69

#t <- 1975 + 1:10
#y <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)


# Simulated data

t <- 1:10
y <- rpois(length(t), exp(7 + 0.2 * t))


# Poisson log-linear model (slightly different from in Gelman)
# Y_i ~ Pois(exp(alpha + beta * t_i)), i=1,...,n


# Compute posterior density (up to proportionality constant) using "uniform prior"

post <- function(a, b) {

  log.prior <- function(a, b) 0
  o <- exp(sum(dpois(y, exp(a + b * t), log=TRUE)) + log.prior(a, b))
  return(o)

}


# Draw contour plot of posterior density

lm.out <- lm(log(y) ~ t)
summary(lm.out)
ci <- confint(lm.out, level=0.99)
A <- seq(ci[1,1], ci[1,2], len=100)
B <- seq(ci[2,1], ci[2,2], len=100)
PP <- 0 * outer(A, B)
for(i in 1:length(A)) {

  for(j in 1:length(B)) PP[i, j] <- post(A[i], B[j])

}
PP <- PP / max(PP)
contour(A, B, PP, nlevels=10, xlab=expression(alpha), ylab=expression(beta), drawlabels=FALSE)


# (Approximate) sampling from the posterior

post.samp <- function(N, A, B, PP) {

  lenA <- length(A); gridA <- A[2] - A[1]
  lenB <- length(B); gridB <- B[2] - B[1]
  L <- lenA * lenB
  I <- sample(L, size=N, replace=TRUE, prob=c(PP))
  aa <- A[(I - 1) %% lenA + 1]
  bb <- B[ceiling(I / lenA)]
  Ajit <- runif(N, -gridA / 2, gridA / 2)
  Bjit <- runif(N, -gridB / 2, gridB / 2)
  out <- matrix(c(aa + Ajit, bb + Bjit), N, 2)
  return(out)

}

AB <- post.samp(1000, A, B, PP)
points(AB[,1], AB[,2], col=2, pch=".", cex=2)


# Estimate the expected number of events in the following year

#est.fit <- exp(AB[,1] + 1986 * AB[,2])
est.fit <- exp(AB[,1] + 11 * AB[,2])
hist(est.fit, breaks=25, freq=FALSE, col="gray", border="white", xlab=expression(exp(alpha + 11~beta)), main="")





