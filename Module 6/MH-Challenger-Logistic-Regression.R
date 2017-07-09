#!/usr/bin/R

# Generic Metropolis-Hastings function
MH <- function(x0, f, dprop, rprop, N, B) {
    ## x0: initializing values for MH algorithm
    ## f: the pdf of posterior; target distribution
    ## dprop: the pdf of the proposal distribution
    ## rprop: the random number gen for proposal
    ## N: Number of samples needed
    ## B: Number of burn in samples
    ## function generates N+B samples then drops B

    x <- matrix(NA, nrow = N + B, ncol = length(x0)) 
    fx <- rep(NA, N + B)
    x[1,] <- x0 ## add the initializing values
    fx[1] <- f(x0) ## initiate the M-H algorithm
    ct <- 0 ## counter for acceptance
    
    for(i in 2:(N + B)){
        u <- rprop(x[i-1,])
        fu <- f(u)
        ## use log values for computational purposes
        r <- log(fu) + log(dprop(x[i-1,])) - log(fx[i-1]) - log(dprop(u))
        R <- min(exp(r), 1)
        if(runif(1) <= R){
            ct <- ct + 1 ## update counter
            x[i,] <- u
            fx[i] <- fu
        } else {
        x[i,] <- x[i-1,]
        fx[i] <- fx[i-1]
        }
    }
    return(list(x=x[-(1:B),], fx=fx[-(1:B)], rate=ct / (N + B)))
}



# Posterior distribution
dpost <- function(theta){
    ## density of Y is binomial/bernoulli
    a <- theta[1]
    b <- theta[2]
    p <- 1 - 1 / (1 + exp(a + b * x)) ## logistic CDF
    lik <- exp(sum(dbinom(y, size=1, prob=p, log=TRUE)))
    dprior <- exp(a) * exp(-exp(a) / b.mme) * 1/ b.mme ## density of prior
    return(lik * dprior)
}


# Proposal distribution (independent proposal, so "theta0" is not used)
dprop <- function(theta){
    ## ignore theta0
    a <- theta[1]
    b <- theta[2]
    pr1 <- exp(a) * exp(-exp(a) / b.mme) / b.mme
    pr2 <- dnorm(b, b.mle, sqrt(var.b.mle))
    return(pr1 * pr2)
}




rprop <- function(theta0){
    ## independent proposals for a and b
    #a <- log(rexp(1, 1 / b.mme))
    a <- log(rexp(1, 1 / b.mme)) ## log for computational purposes
    b <- rnorm(1, b.mle, sqrt(var.b.mle))
    return(c(a, b))
}





# O-ring failure data
## Y = 1 (failure)
y <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0)
x <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70, 70, 70, 72, 73, 75, 75,
76, 76, 78, 79, 81)




# Preliminary output from ML estimation
fit <- glm(y ~ x, family=binomial(logit))
print(summary(fit))
a.mle <- as.numeric(fit$coefficients[1])
b.mle <- as.numeric(fit$coefficients[2])
var.a.mle <- summary(fit)$cov.scaled[1, 1]
var.b.mle <- summary(fit)$cov.scaled[2, 2]
b.mme <- exp(a.mle + 0.577216)






# Run Metropolis-Hastings
N <- 10^5
B <- 10^4
mh.out <- MH(c(a.mle, b.mle), dpost, dprop, rprop, N, B)
alpha.mh <- mh.out$x[,1]
beta.mh <- mh.out$x[,2]

print(summary(mh.out$x))

par(mfrow = c(1,2))
hist(alpha.mh, freq=FALSE, col="gray", border="white", xlab=expression(alpha))
plot(alpha.mh, type="l", col="gray", xlab="Iteration", ylab=expression(alpha))
lines(1:N, cumsum(alpha.mh) / (1:N))

hist(beta.mh, freq=FALSE, col="gray", border="white", xlab=expression(beta))
plot(beta.mh, type="l", col="gray", xlab="Iteration", ylab=expression(beta))
lines(1:N, cumsum(beta.mh)/(1:N))

p65 <- 1 - 1 / (1 + exp(alpha.mh + beta.mh * 65))
p45 <- 1 - 1 / (1 + exp(alpha.mh + beta.mh * 45))

hist(p65, freq=FALSE, col="gray", border="white", xlab="p(65)", main="")
hist(p45, freq=FALSE, col="gray", border="white", xlab="p(45)", main="")




