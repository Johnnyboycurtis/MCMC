
##--------------------------------------------------------------------##
## Stat 591                                                           ##
## Author: Ryan Martin (rgmartin@math.uic.edu)                        ##
## Date: 08/30/2013                                                   ##
## R code for likelihood & bootstrap inference on the gamma variance  ##
##--------------------------------------------------------------------##


## Generic function for the Newton (or Newton-Raphson) optimization method

# Note: The optimization methods built in to R (e.g., 'optim') are better 
# than this implementation of the Newton method given here, which is just 
# for illustration purposes (so students can see how to write their own
# implementation of Newton's method).  With this version, you might find 
# that if the initial guess is poor, then the method will not converge to 
# a solution.  The built-in methods in R are faster and more stable.

newton <- function(f, df, x0, eps=1e-08, maxiter=1000, ...) {

  if(!exists("ginv")) library(MASS)
  x <- x0
  t <- 0

  repeat {

    t <- t + 1
    x.new <- x - as.numeric(ginv(df(x, ...)) %*% f(x, ...))
    if(mean(abs(x.new - x)) < eps | t >= maxiter) {

      if(t >= maxiter) warning("Maximum number of iterations reached!")
      break

    }
    x <- x.new

  }

  out <- list(solution=x.new, value=f(x.new, ...), iter=t)
  return(out)

}



## Auxiliary functions

# Derivative of the gamma log-likelihood

dl <- function(theta, X, n) {

  alpha <- theta[1]
  beta <- theta[2]
  o1 <- -n * log(beta) - n * digamma(alpha) + sum(log(X))
  o2 <- -n * alpha / beta + n * mean(X) / beta**2
  return(c(o1, o2))

}

# Second derivative of the gamma log-likelihood

ddl <- function(theta, X, n) {

  alpha <- theta[1]
  beta <- theta[2]
  o11 <- -n * trigamma(alpha)
  o12 <- -n / beta
  o22 <- -n * (2 * mean(X) / beta**3 - alpha / beta**2)
  return(matrix(c(o11, o12, o12, o22), 2, 2, byrow=TRUE))

}

# Gamma method of moments estimation

gamma.mme <- function(X) {

  m <- mean(X)
  v <- var(X)
  return(c(m**2 / v, v / m))

}


## Gamma maximum likelihood estimation

X <- c(6.00, 5.98, 7.81, 6.77, 10.64, 13.63, 10.53, 8.92, 3.77, 5.78,
       6.32, 4.44, 5.34, 1.54, 4.42, 4.71, 6.12, 2.57, 9.47, 9.03)
n <- length(X)
mme <- gamma.mme(X)                             # method of moments estimator
mle <- newton(dl, ddl, mme, X=X, n=n)$solution  # maximum likelihood estimator
print(rbind(mme=mme, mle=mle))


## Asymptotic likelihood inference on the variance (= alpha * beta**2)

J <- -ddl(mle, X, n)                                                # observed information
gr <- matrix(c(mle[2]**2, 2 * mle[1] * mle[2]), ncol=1)             # gradient of variance
dt.var <- as.numeric(t(gr) %*% ginv(J) %*% gr)                      # delta theorem variance
mle.ci <- mle[1] * mle[2]**2 + qnorm(c(0.05, 0.95)) * sqrt(dt.var)  # 90% CI for variance
print(mle.ci)


## Bootstrap inference on the variance

B <- 2500


# Nonparametric bootstrap

var.nboot <- numeric(B)
for(b in 1:B) var.nboot[b] <- var(sample(X, size=n, replace=TRUE))
nboot.ci <- as.numeric(quantile(var.nboot, c(0.05, 0.95)))  # 90% (npar) bootstrap CI for variance
print(nboot.ci)


# Parametric bootstrap

var.pboot <- numeric(B)
for(b in 1:B) {

  XX <- rgamma(n, shape=mle[1], scale=mle[2])
  mmeb <- gamma.mme(XX)
  mleb <- newton(dl, ddl, mmeb, X=XX, n=n)$solution
  var.pboot[b] <- mleb[1] * mleb[2]**2

}
pboot.ci <- as.numeric(quantile(var.pboot, c(0.05, 0.95)))  # 90% (par) bootstrap CI for variance
print(pboot.ci)


## Comparison of "sampling distributions"

hist(var.nboot, freq=FALSE, col="gray", border="white", xlab="Variance", xlim=c(0,25), ylim=c(0,0.15), main="")
hist(var.pboot, freq=FALSE, add=TRUE)
curve(dnorm(x, mle[1] * mle[2]**2, sqrt(dt.var)), add=TRUE, lwd=2)
legend(x="topright", inset=0.05, fill=c("gray", "white"), c("npar", "par"))




