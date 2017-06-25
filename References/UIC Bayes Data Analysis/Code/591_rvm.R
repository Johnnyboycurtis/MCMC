
##---------------------------------------------------------------------------##
## Stat 591                                                                  ##
## Author: Ryan Martin (rgmartin@uic.edu)                                    ##
## R code for accept-reject sampling from a von Mises distribution           ##
## See Lange (2010), Chapter 22.5, in particular, Exercise 22.17             ##
##---------------------------------------------------------------------------##

# A bisection routine to be called in rvm()

bisection <- function (f, a, b, eps=1e-08, maxiter=1000, verbose=FALSE, ...) {

  x <- (a + b) / 2
  t <- 0
  if(verbose) cat("\n")
  
  repeat {

    t <- t + 1
    if(f(a, ...) * f(x, ...) <= 0) b <- x else a <- x
    x.new <- (a + b) / 2
    if(verbose) cat("[", t, "]", x.new, f(x.new, ...), "\n")
    if(abs(x.new - x) < eps | t >= maxiter) {

      if(t >= maxiter) warning("Maximum number of iterations reached!")
      break

    }
    x <- x.new

  }
  if(verbose) cat("\n")
  out <- list(solution=x.new, value=f(x.new, ...), iter=t)
  return(out)

}

# von-Mises distribution density function (normalized and unnormalized)

dvm <- function(x, kappa=1) {

  if(kappa < 0) stop("kappa must be non-negative!")
  d <- exp(kappa * cos(x)) * ifelse(x > 0 & x < 2 * pi, 1, 0) / (2 * pi * besselI(kappa, 0))
  return(d)

}

dvm0 <- function(x, kappa=1) {

  if(kappa < 0) stop("kappa must be non-negative!")
  d <- exp(kappa * cos(x)) * ifelse(x > 0 & x < 2 * pi, 1, 0) # / (2 * pi * besselI(kappa, 0))
  return(d)

}


# Visualize the distribution to be sampled from

curve(dvm0(x, kappa=2), xlim=c(0.1, 2 * pi - 0.1))


# Function to sample from von-Mises distribution using acceptance-rejection method.
# The goal is to use a pair of opposite oriented exponential densities as majorants.
# Most of the code is devoted to choosing the "best" exponential majorants.  

rvm <- function(n, kappa=1, pic=FALSE) {

  if(kappa < 0) stop("kappa must be non-negative!")
  g <- function(x) kappa * x * sin(x) - 1
  x.left <- bisection(f=g, a=0, b=pi)$solution
  x.right <- 2 * pi - x.left
  lambda <- kappa * sin(x.left)
  k <- exp(kappa * (x.left * cos(x.left) + sin(x.left))) / (2 * pi * kappa * sin(x.left) * besselI(kappa, 0))
  maj <- function(x) {

    out <- k * lambda * exp(-lambda * ifelse(x < pi, x, 2 * pi - x))
    return(out)

  }
  if(pic) {

    curve(dvm(x, kappa), xlim=c(0.1, 2*pi-0.1), lwd=2, col="red")
    curve(maj(x), add=TRUE, lwd=2, col="blue")
    legend(x="top", inset=0.05, lwd=2, col=c("red", "blue"), c("True", "Majorant"))

  }
  r <- function(x) dvm0(x, kappa) / maj(x)
  nsamp <- 0
  count <- 0
  X <- rep(NA, n)
  while(nsamp < n) {

    count <- count + 1
    x <- rexp(1, rate=lambda)
    if(runif(1) <= r(x)) {

      nsamp <- nsamp + 1
      X[nsamp] <- ifelse(runif(1) <= 0.5, x, 2 * pi - x)

    }

  }
  AcceptRate <- n / count
  print(cbind(x.left=x.left, lambda=lambda, k=k, AcceptRate=AcceptRate))
  return(X)

}


# Example

vm.out <- rvm(n=1000, kappa=2, pic=TRUE)
hist(vm.out, breaks=20, freq=FALSE, xlab="x", col="gray", border="white", main="")
curve(dvm(x, kappa=2), xlim=c(0.1, 2 * pi - 0.1), add=TRUE, lwd=2)





