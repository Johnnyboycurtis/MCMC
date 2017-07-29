library(coda)

# This function takes the prior mode and a percentile and returns the corresponding hyperparameter values for the beta prior
find.beta <- function(mode, percentile, p=0.95, start.b=1, end.b=1000, increment=0.01){
  b <- seq(from=start.b, to=end.b, by=increment)
  a <- (1 + mode * (b - 2)) / (1 - mode)
  i <- which.max(qbeta(p, a, b) < percentile)
  return(list(a=a[i], b=b[i])) 
}
# This function takes the prior mean and a percentile and returns the corresponding hyperparameter values for the normal prior
find.normal <- function(prior.mean, percentile, p=0.95){
  prior.sd <- (qnorm(p) / (percentile - prior.mean))^{-1}
  return(list(mean=prior.mean, sd=prior.sd, precision=1/prior.sd^2)) 
}

# This function takes the prior mode and a percentile of a random variable with a Gamma distribution and returns the corresponding parameters
find.gamma <- function(prior.mode, percentile, p=0.95, start.b=0.001, end.b=1000, increment=0.01){
  b <- seq(from=start.b, to=end.b, by=increment)
  a <- 1 + b * prior.mode
  i <- which.max(qgamma(p, a, b) < percentile)
  return(list(a=a[i], b=b[i])) 
}

# This function takes the prior mode and a percentile of sigma and returns the corresponding parameters of the Gamma distribution for tau 
find.tau.gamma <- function(prior.sigma.mode, sigma.percentile, p=0.95, start.a=1, end.a=1000, increment=0.01){
  a <- seq(from=start.a, to=end.a,  by=increment)
  b <- prior.sigma.mode^2 * (a + 1)
  i <- which.max(qgamma(1-p, a, b) > 1/sigma.percentile^2)
  return(list(a=a[i], b=b[i])) 
}

# Given the value of the mean, this function converts normal percentiles to standard deviations. 
normal.percentile.to.sd <- function(mean.value, percentile, p=0.95){
  return(list(sd=(percentile-mean.value)/ qnorm(p)))  
}


# This is function provides a simple summary of a unimodal posterior distribution
summary.Bayes <- function(posterior.sample, HPD.p=0.95){
  basic.stats <- c(summary(posterior.sample), sd(posterior.sample), HPDinterval(as.mcmc(posterior.sample)), HPD.p)
  names(basic.stats) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max", "SD", "HPD-Lower", "HPD-Upper", "Level")
  return(basic.stats)
}

library(MASS)

# Gibbs sampler for linear regression model with independent proper priors. 
Gibbs.lm <- function(net.iterations, n.thin=1, burn.in=10000, start.beta,  X, y, a, b, prior.mean.beta, prior.cov.beta){
  N <- length(y)
  r <- ncol(X)
  prior.cov.beta.inv <- solve(prior.cov.beta) 
  
  # Terms to be computed once and stored
  beta.prior.term <- prior.cov.beta.inv %*% prior.mean.beta 
  beta.hat <- solve(t(X)%*%X) %*% t(X) %*% y
  e <- y - X %*% beta.hat 
  SSE <- t(e) %*% e 
  a.pos <- a + N/2.0
  tXX <- t(X) %*% X
  tXy <- t(X) %*% y 
  
  n.iterations <- burn.in + net.iterations 
  beta <- matrix(NA, nrow = n.iterations, ncol = r)
  tau <- rep(NA, length = n.iterations)
  temp_beta <- start.beta
  V_inv <- matrix(NA, nrow = r, ncol = r)
  
  for(i in 1:n.iterations){
    for(j in 1:n.thin){ 
      temp_tau <- rgamma(1, shape = a.pos, rate = b + (SSE + t(temp_beta - beta.hat) %*% tXX %*% (temp_beta - beta.hat)) / 2.0)
      V_inv <- temp_tau * tXX + prior.cov.beta.inv
      V <- chol2inv(chol(V_inv))
      temp_beta <- V %*% (temp_tau * tXy + beta.prior.term) + t(chol(V)) %*% rnorm(r)
      # temp_beta <- mvrnorm(n = 1, mu= V %*% (temp_tau * tXy + beta.prior.term),  Sigma= V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    }
    beta[i , ] <- temp_beta
    tau[i] <- temp_tau
  }
  sigma <- 1 / sqrt(tau)
  beta <- beta[(burn.in+1):n.iterations, ]
  sigma <- sigma[(burn.in+1):n.iterations]
  return( list(beta=beta, sigma=sigma))  
}

# This function takes the MCMC output and a vector of covariates, and returns the plot of the predictive density. 
Gibbs.predict <- function(MCMC, x.f, lognormal=FALSE, main="", xlim){
  y.f.sample <- MCMC$beta %*% x.f + MCMC$sigma * rnorm(nrow(MCMC$beta)) # This replaced a `for' loop.
  if(lognormal==TRUE){
    y.f.sample <- exp(y.f.sample)
  }
  summary.Bayes(as.vector(y.f.sample))
  # predictive density of future salary for someone with x.f covariate values 
  plot(density(y.f.sample), type="l", main=main, xlab="new observation", ylab="density", xlim=xlim)
  lines(HPDinterval(as.mcmc(y.f.sample)), c(0,0), col="red", cex=2, lty=3)  # add PI line
  legend("topright", legend="95% HPD interval", lty=2, col="red", cex=0.5)
}


# Gibbs sampler for linear mixed regression model with improper priors.
GibbslmeImproper<- function(net.iterations, n.thin=1, burn.in, start.theta,  X, Z, y, a = c(0, -1/2), b = c(0,0)){
  N <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  W <- cbind(X, Z)
  a_pos <- c( a[1] + N/2.0, a[2] + q/2.0)
  
  tXX <- t(X) %*% X
  tXZ <- t(X) %*% Z
  tZX <- t(Z) %*% X
  tZZ <- t(Z) %*% Z
  tXy <- t(X) %*% y 
  tZy <- t(Z) %*% y
  
  
  n.iterations <- burn.in + net.iterations 
  thetas <- matrix(0, nrow = n.iterations, ncol = {p+q})
  lambdas <- matrix(0, nrow = n.iterations, ncol = 2)
  temp_thetas <- start.theta
  temp_lambdas <- c(0,0)
  I.q <- diag(q)
  eta <- rep(0, p+q)
  V_inv <- matrix(0 , nrow = p + q, ncol = p + q)
  
  for( i in 1:n.iterations){
    for(j in 1:n.thin){
      diff = y - W %*% temp_thetas
      temp_lambdas[1] <- rgamma(1, shape = a_pos[1], rate = b[1] + t(diff) %*% diff / 2.0)
      diffu = temp_thetas[{p+1}:{p+q}]
      temp_lambdas[2] <- rgamma(1, shape = a_pos[2], rate = b[2] + t(diffu) %*% diffu / 2.0)
      
      V_inv[1:p, 1:p] <- temp_lambdas[1] * tXX
      V_inv[1:p, {p+1}:{p+q}] <-  temp_lambdas[1] * tXZ
      V_inv[{p+1}:{p+q}, 1:p] <-  temp_lambdas[1] * tZX
      V_inv[{p+1}:{p+q}, {p+1}:{p+q}] <-  temp_lambdas[1] * tZZ + temp_lambdas[2] * I.q
      V <- chol2inv(chol(V_inv))
      eta[1:p] <- temp_lambdas[1] * tXy 
      eta[{p+1}:{p+q}] <- temp_lambdas[1] * tZy
      temp_thetas <- V %*% eta + t(chol(V)) %*% rnorm(p+q)
    }
    thetas[i , ] <- temp_thetas
    lambdas[i , ] <- temp_lambdas
  }
  sigmas <- 1 / sqrt(lambdas)
  thetas <- thetas[{burn.in+1}:n.iterations, ]
  sigmas <- sigmas[{burn.in+1}:n.iterations, ]
  return( list( beta = thetas[ , 1:p], group = thetas[ , {p+1}:{p+q}], sigma = sigmas) )  
}





