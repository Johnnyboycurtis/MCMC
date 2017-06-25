
##--------------------------------------------------------------------##
## Stat 591                                                           ##
## Author: Ryan Martin (rgmartin@uic.edu)                             ##
## Date: 10/18/2013                                                   ##
## R code for Bayes inference in the many-normal-mean problem         ##
##--------------------------------------------------------------------##


# Model: X_i ~ N(mu_i, 1), independent, i=1,...,n.
# Goal is inference on psi = sqrt(sum(mu**2)), length of the mu vector.
# (In class I worked with eta = psi**2, but there's no difference.)


# Plot posterior density of psi, with credible interval

mnm.post.plot <- function(n, mu, min, max, alpha=0.05) {

  if(length(mu) > 1 && length(mu) != n) stop("mu is not of the right dimension!")
  if(length(mu) == 1) mu <- rep(mu, n)
  psi <- sqrt(sum(mu * mu))
  X <- rnorm(n) + mu
  ssX <- sum(X * X)
  P <- seq(min, max, length=300)
  dpost <- 2 * P * dchisq(P**2, df=n, ncp=ssX) 
  plot(0, 0, type="n", xlim=range(P), ylim=range(dpost), xlab=expression(psi), ylab="Posterior")
  abline(v=psi, col="gray", lwd=2)
  abline(h=0)
  qtile <- c(alpha / 2, 1 - alpha / 2)
  ba.int <- sqrt(qchisq(qtile, df=n, ncp=ssX))
  points(x=ba.int, y=c(0, 0), pch=c("(", ")"))
  lines(P, dpost, lwd=2)
  return(ssX)

}


# Simulate data and evaluate coverage probability of credible intervals

mnm.bayes.sim <- function(n, mu, reps=1e4, alpha=0.05) {

  if(length(mu) > 1 && length(mu) != n) stop("mu is not of the right dimension!")
  if(length(mu) == 1) mu <- rep(mu, n)
  psi <- sqrt(sum(mu * mu))
  cvg <- 0
  for(r in 1:reps) {

    X <- mu + rnorm(n)
    ssX <- sum(X * X)
    b <- pchisq(psi**2, df=n, ncp=ssX) 
    cvg <- cvg + (b > alpha / 2 && b < 1 - alpha / 2) / reps

  }
  return(cvg)

}


# TO RUN:
# mnm.post.plot(n=10, mu=1, 0, 20)
# mnm.sim(n=10, mu=1)


