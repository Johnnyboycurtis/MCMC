#!/usr/bin/R

# Bayesian Analysis
## O-ring Challenger Data { .selectable }

  

failure <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 
             0, 0, 0, 1, 0, 0, 0, 0, 0)
temperature <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 
                 69, 70, 70, 70, 70, 72, 73, 75, 75,
                 76, 76, 78, 79, 81)

df = data.frame(failure, temperature)
head(df)


## O-ring Challenger Data { .selectable }
#The frequentist logistic regression

fit = glm(formula = failure ~ temperature, data = df, 
          family = binomial(link = "logit"))

summary(fit)



## O-ring Challenger Data { .selectable }
## Output from ML estimation from Logistic Regression
## MLEs make great starting valules
a.mle <- as.numeric(fit$coefficients[1])
b.mle <- as.numeric(fit$coefficients[2])
var.a.mle <- summary(fit)$cov.scaled[1, 1]
var.b.mle <- summary(fit)$cov.scaled[2, 2]

b.hyper <- exp(a.mle + 0.577216) ## hyper parameter
## 0.577216 is "Euler’s constant"


## Let's set up the posterior and proposal distributions


## setting up functions

# Posterior distribution
dPosterior <- function(theta, y = failure, x = temperature){
## density of Y is binomial/bernoulli
a <- theta[1]
b <- theta[2]
p <- 1 - 1 / (1 + exp(a + b * x)) ## logistic CDF
lik <- exp(sum(dbinom(y, size=1, prob=p, log=TRUE)))
dprior <- exp(a) * exp(-exp(a) / b.hyper) * 1/b.hyper ## density of prior
return(lik * dprior)
}




## O-ring Challenger Data { .selectable }


# Proposal distribution (independent proposal, so "theta0" is not used)
dProposal <- function(theta){
## ignore theta0
a <- theta[1]
b <- theta[2]
# a <- log(rexp(1, 1 / b.hyper)) ## remember, log for computational purposes
# try: exp(rexp(1, 1/b.hyper)) ## Inf
pr1 <- exp(a) * exp(-exp(a) / b.hyper) * 1/b.hyper
pr2 <- dnorm(b, b.mle, sqrt(var.b.mle))
return(pr1 * pr2)
}



## O-ring Challenger Data { .selectable }


rProposal <- function(theta0){
## independent proposals for a and b
#a <- log(rexp(1, 1 / b.hyper))
a <- log(rexp(1, 1 / b.hyper)) ## log for computational purposes
b <- rnorm(1, b.mle, sqrt(var.b.mle))
return(c(a, b))
}



## Metropolis-Hastings set up to run
# Run Metropolis-Hastings
N = 1000000
BurnIn =  5000




start <- Sys.time()

## x0: initializing values for MH algorithm
## f: the pdf of posterior; target distribution
## dProposal: the pdf of the proposal distribution
## rProposal: the random number gen for proposal
## N: Number of samples needed
## BurnIn: Number of burn in samples
## function generates N+B samples then drops B
x0 <- c(a.mle, b.mle)
x <- matrix(NA, nrow = N + BurnIn, ncol = length(x0)) ## empty matrix
fx <- rep(NA, N + BurnIn)
x[1,] <- x0 ## add the initializing values
fx[1] <- dPosterior(x0) ## initiate the M-H algorithm
accept <- 0 ## counter for acceptance

for(i in 2:(N + BurnIn)){
u <- rProposal(x[i-1,])
fu <- dPosterior(u)
## use log values for computational purposes
r <- log(fu) + log(dProposal(x[i-1,])) - log(fx[i-1]) - log(dProposal(u))
Rho <- min(exp(r), 1)
if(runif(1) <= Rho){
accept <- accept + 1 ## update counter
x[i,] <- u
fx[i] <- fu
} else {
x[i,] <- x[i-1,]
fx[i] <- fx[i-1]
}
}



print("acceptance rate: ")
print(accept/(N+BurnIn))

MH.Results <- x[-(1:BurnIn), ]
alphaMH <- MH.Results[,1]
betaMH <- MH.Results[,2]


summary(MH.Results) ## summary of parameters



end <- Sys.time()
print(paste("start: ", start))
print(paste("end: ", end))
print(paste("Total Time to run MH:", (end-start)))




par(mfrow = c(1,2))
hist(alphaMH, freq=FALSE, col="gray", border="white", xlab=expression(alpha), main = "")
plot(alphaMH, type="l", col="gray", xlab="Iteration", ylab=expression(alpha), main = "")
lines(1:N, cumsum(alphaMH) / (1:N))


par(mfrow = c(1,2))
hist(betaMH, freq=FALSE, col="gray", border="white", xlab=expression(beta), main = "")
plot(betaMH, type="l", col="gray", xlab="Iteration", ylab=expression(beta), main = "")
lines(1:N, cumsum(betaMH)/(1:N))



