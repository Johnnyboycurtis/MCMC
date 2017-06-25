
##---------------------------------------------------------------------------##
## Stat 591                                                                  ##
## Author: Ryan Martin (rgmartin@uic.edu)                                    ##
## R code implementation of the Metropolis-Hastings algorithm                ##
## Two illustrations                                                         ##
##---------------------------------------------------------------------------##


## General code for the Metropolis-Hastings algorithm

# INPUT:
#    x0 = starting point of the chain (possibly a vector)
#     f = target density function
# dprop = density function for the proposal distribution
# rprop = sampling from the proposal distribution
#     N = Monte Carlo sample size
#     B = burn-in length

# OUTPUT:
#    x = matrix of MH samples (one in each row)
#   fx = values of the posterior density at the sampled points
# rate = proportion of times the proposal draw was accepted

mh <- function(x0, f, dprop, rprop, N, B) {

  x <- matrix(NA, N + B, length(x0))
  fx <- rep(NA, N + B)
  x[1,] <- x0
  fx[1] <- f(x0)
  ct <- 0
  for(i in 2:(N + B)) {

    u <- rprop(x[i-1,])
    fu <- f(u)
    r <- log(fu) + log(dprop(x[i-1,], u)) - log(fx[i-1]) - log(dprop(u, x[i-1,]))
    R <- min(exp(r), 1)

    if(runif(1) <= R) {

      ct <- ct + 1
      x[i,] <- u
      fx[i] <- fu

    } else {

      x[i,] <- x[i-1,]
      fx[i] <- fx[i-1]

    }

  }
  out <- list(x=x[-(1:B),], fx=fx[-(1:B)], rate=ct / (N + B))
  return(out)

}


## Example: Simulate from standard (independent) bivariate normal (Gelman et al, 2004, p. 290)

# Density function for proposal distribution

dprop <- function(x, theta) {

  scale <- 0.2
  out <- exp(sum(dnorm(as.numeric(x - theta), 0, scale, log=TRUE)))
  return(out)

}

# Sampling function for the proposal distribution

rprop <- function(theta) {

  scale <- 0.2
  theta.new <- as.numeric(theta) + scale * rnorm(length(theta))
  return(theta.new)

}


# Target posterior density

f <- function(theta) return(exp(sum(dnorm(as.numeric(theta), log=TRUE))))


# Reproduce results in Gelman et al, 2004, Figure (c), page 286

theta1 <- seq(-3, 3, len=100)
theta2 <- seq(-3, 3, len=100)
dtheta <- 0 * outer(theta1, theta2)

plot(x=0, y=0, xlim=c(-3, 3), ylim=c(-3, 3), type="n", xlab=expression(theta[1]), ylab=expression(theta[2]))
for(i in 1:4) {

  out <- mh(runif(2, -1, 1), f, dprop, rprop, 1000, 1000)
  points(out$x[,1], out$x[,2], pch=".", cex=2, col=i)

}


## Example: Simulate from posterior in a trigonometric location problem (Givens & Hoeting)

X <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96, 2.53,
       3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52, 2.50)
       
lik <- function(theta) exp( sum( log( 1 - cos(X - theta) ) ) )
a <- 0.5
dprop <- function(theta, theta0) dunif(theta, theta0-a, theta0+a)
rprop <- function(theta0) runif(1, theta0-a, theta0+a)
N <- 10000
B <- 1000
theta.mcmc <- mh(runif(1, -pi, pi), lik, dprop, rprop, N, B)
hist(theta.mcmc$x, breaks=25, freq=FALSE, col="gray", border="white", xlab=expression(theta), main="Posterior Sample")
plot(theta.mcmc$x, type="l", col="gray", xlab="Iteration (t)", ylab=expression(theta[t]))
lines(1:N, cumsum(theta.mcmc$x)/(1:N))
print(mean(theta.mcmc$x))
print(quantile(theta.mcmc$x, c(0.05, 0.95)))

