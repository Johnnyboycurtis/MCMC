
##--------------------------------------------------------------------##
## Stat 591                                                           ##
## Author: Ryan Martin (rgmartin@math.uic.edu)                        ##
## Date: 09/11/2013                                                   ##
## R code for Bayesian inference in binomial-beta model               ##
##--------------------------------------------------------------------##


## Bayes model:
## X | theta ~ Bin(n, theta)
## theta ~ Beta(a, b) [conjugate prior]


## Plots of the prior and posterior

op <- par(mfrow=c(3,3), mar=c(4.2, 4.2, 3, 1))
n <- c(10, 25, 50)
hpar <- matrix(c(1, 1, 0.5, 0.5, 10, 25), nrow=3, byrow=TRUE)
mle <- 0.6
for(i in 1:3) {
	
  a <- hpar[i,1]
  b <- hpar[i,2]
  a0 <- a + n[3] * mle
  b0 <- b + n[3] * (1 - mle)
  top <- dbeta((a0 - 1) / (a0 + b0 - 2), a0, b0)
  for(j in 1:3) {
  	
  	ax <- a + n[j] * mle
  	bx <- b + n[j] * (1 - mle)
  	mn <- paste("a=", a, " b=", b, " n=", n[j], sep="")
  	plot(0, 0, type="n", ylab="Posterior pdf", xlab=expression(theta), ylim=c(0, top), xlim=c(0,1), main=mn)
  	curve(dbeta(x, shape1=a, shape2=b), ylim=c(0, top), col="gray", add=TRUE)
  	curve(dbeta(x, shape1=ax, shape2=bx), add=TRUE)
  	
  }
	
}
par(op)


## Credible intervals for theta

# Auxiliary function: the bisection method for solving equations

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



# Equi-tail 90% credible interval

cred.int <- function(x, n, a, b, alpha=0.10) {
	
  ax <- a + x
  bx <- b + n - x
  out <- numeric(2)
  cut <- c(alpha / 2, 1 - alpha / 2)
  for(i in 1:2) {
  	
    F <- function(u) pbeta(u, ax, bx) - cut[i]
    out[i] <- bisection(F, 0, 1)$solution
  
  }
  return(out)
	
}


# Highest posterior density 90% credible interval

cred.int.hpd <- function(x, n, a, b, alpha=0.10) {
	
  ax <- a + x
  bx <- b + n - x
  mode <- (ax - 1) / (ax + bx - 2)
  out <- numeric(2)
  f <- function(r) {
  	
  	o <- numeric(2)
  	g <- function(u) dbeta(u, ax, bx) - r
  	o[1] <- bisection(g, 0, mode)$solution
  	o[2] <- bisection(g, mode, 1)$solution
  	return(list(low=o[1], up=o[2], pr=pbeta(o[2], ax, bx) - pbeta(o[1], ax, bx)))
  	
  }
  ff <- function(r) f(r)$pr
  F <- function(r) ff(r) - (1 - alpha)
  R <- bisection(F, 1e-3, dbeta(mode, ax, bx) - 1e-3)$solution
  fR <- f(R)
  out[1] <- fR$low
  out[2] <- fR$up
  return(out)
	
}


# Example:

cred.int(x=7, n=12, a=1, b=1)
cred.int.hpd(x=7, n=12, a=1, b=1)


# Simulation: frequentist coverage probability check 

cred.int.sim <- function(theta, reps, n, a, b, alpha=0.10) {
	
  cp <- 0
  for(r in 1:reps) {
  	
  	X <- rbinom(1, size=n, prob=theta)
  	ci <- cred.int(X, n, a, b, alpha)   # can replace 'cred.int' with 'cred.int.hpd'
  	cp <- cp + (ci[1] <= theta && ci[2] >= theta) 
  	
  }
  cp <- cp / reps
  return(cp)
	
}

# TO RUN:
cred.int.sim(0.25, 1000, 10, 1, 1)


