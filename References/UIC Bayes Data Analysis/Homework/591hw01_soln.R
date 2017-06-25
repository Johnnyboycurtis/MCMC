
##--------------------------------------------------------------------##
## Stat 591                                                           ##
## Author: Ryan Martin (rgmartin@uic.edu)                             ##
## Date: 09/11/2013                                                   ##
## R code solutions for Homework 01                                   ##
##--------------------------------------------------------------------##

## Kappa 

# Function kappa = g(theta)

g <- function(u) {
	
  o <- (u[1] + u[2]) * (u[1] + u[3]) + (u[3] + u[4]) * (u[2] + u[4])
  fn <- ((u[1] + u[4] - o) / (1 - o))
  return(fn)
	
}

# (Numerical) gradient of the kappa function

dg <- function(u) {
	
  eps <- 1e-08
  gr <- numeric(4)
  for(i in 1:4) {
  	
  	uu <- u
  	uu[i] <- u[i] + eps
  	gr[i] <- (g(uu) - g(u)) / eps
  	
  }	
  return(matrix(gr, nrow=4))
	
}


## Covariance matrix of the MLE (= multinomial covariance / n**2)

cmat <- function(u, n) {
	
  o <- diag(u) - outer(u, u)
  return(o / n)
	
}


## Function for bootstrap distribution of kappa.mle

boot <- function(X, n, B=3000) {
	
  kap <- numeric(B)
  mle <- X / n
  for(b in 1:B) kap[b] <- g(as.numeric(table(sample(1:4, size=n, prob=mle, replace=TRUE))) / n)
  return(kap)
	
}


# TO RUN:

X <- c(22, 15, 33, 30)
n <- sum(X)
mle <- X / n
kappa.mle <- g(mle)
dg.mle <- dg(mle)
kappa.var <- t(dg.mle) %*% cmat(mle, n) %*% dg.mle
kappa.mle.ci <- kappa.mle + c(-1, 1) * qnorm(0.90) * sqrt(kappa.var)
kappa.boot <- boot(X, n)
kappa.boot.ci <- as.numeric(quantile(kappa.boot, c(0.05, 0.95)))

print(kappa.var)
print(kappa.mle.ci)
print(kappa.boot.ci)

hist(kappa.boot, freq=FALSE, col="gray", border="white", main="", xlab=expression(hat(kappa)))
curve(dnorm(x, kappa.mle, sqrt(kappa.var)), add=TRUE)




