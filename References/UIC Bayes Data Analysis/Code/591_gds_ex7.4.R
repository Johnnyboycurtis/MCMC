
##---------------------------------------------------------------------------##
## Stat 591                                                                  ##
## Author: Ryan Martin (rgmartin@uic.edu)                                    ##
## R code for Monte Carlo in Example 7.4 in [GDS]                            ##
##---------------------------------------------------------------------------##

# Model: X1,...,Xn ~ N(mu, sigma2), sigma2 known
# Prior: mu ~ Cauchy(loc, scale)


# Data generation

loc <- 0
scale <- 1
mu <- 5
sigma <- 2
n <- 5
Xbar <- rnorm(1, mu, sigma / sqrt(n))


# Exact value of posterior mean via numerical integration

num <- function(u) u * dnorm(u, Xbar, sigma / sqrt(n)) * dt((u - loc) / scale, df=1)
den <- function(u) dnorm(u, Xbar, sigma / sqrt(n)) * dt((u - loc) / scale, df=1)
exact <- integrate(num, -Inf, Inf)$value / integrate(den, -Inf, Inf)$value


# Monte Carlo approach 1: Simulate from normal

mc1 <- function(M) {

  tt <- rnorm(M, Xbar, sigma / sqrt(n))
  mc.num <- sum(tt * dt((tt - loc) / scale, df=1))
  mc.den <- sum(dt((tt - loc) / scale, df=1))
  return(mc.num / mc.den)

}


# Monte Carlo approach 2: Simulate from Cauchy

mc2 <- function(M) {

  tt <- loc + scale * rt(M, df=1)
  mc.num <- sum(tt * dnorm(tt, Xbar, sigma / sqrt(n)))
  mc.den <- sum(dnorm(tt, Xbar, sigma / sqrt(n)))
  return(mc.num / mc.den)

}


# Compare results

M <- c(500, 1000, 5000, 10000)
cbind(exact=exact, MC1=sapply(M, mc1), MC2=sapply(M, mc2))

