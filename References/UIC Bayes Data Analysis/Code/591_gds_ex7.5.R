
##---------------------------------------------------------------------------##
## Stat 591                                                                  ##
## Author: Ryan Martin (rgmartin@uic.edu)                                    ##
## R code for importance sampling in Example 7.5 in [GDS]                    ##
##---------------------------------------------------------------------------##

# Model: X1,...,Xn ~ N(mu, sigma2), both mu and sigma2 unknown
# Prior: mu ~ DoubleExp and sigma2 ~ HalfStudent-t


# Data generation

mu <- 5
sigma <- 2
n <- 5
X <- rnorm(n, mu, sigma)
Xbar <- mean(X)
S <- (1 - 1 / n) * var(X)


# Likelihood and prior

L <- function(mu, sigma2) exp(sum(dnorm(X, mu, sqrt(sigma2), log=TRUE)))
dprior <- function(mu, sigma2) exp(-abs(mu)) / (1 + sigma2)**2


# Importance function and weight

g <- function(mu, sigma2) {
	
  v <- sigma2
  D <- (n / 2) * ((mu - Xbar)**2 + S)
  g2 <- dt(mu, df=n+1)
  g1 <- (1 / v)**(n / 2) * exp(- 1 / D / v)
  return(g1 * g2)
	
}

w <- function(mu, sigma2) L(mu, sigma2) * dprior(mu, sigma2) / g(mu, sigma2)


# Importance sampling to approximate posterior mean of mu

M <- 5000
num <- 0
den <- 0
for(m in 1:M) {
	
  mu <- rt(1, n+2)
  D <- (n / 2) * ((mu - Xbar)**2 + S)
  sigma2 <- 1 / rgamma(1, shape=n/2, scale=D)
  ww <- w(mu, sigma2)
  num <- num + mu * ww
  den <- den + ww
	
}
postmean.is <- num / den
print(postmean.is)

