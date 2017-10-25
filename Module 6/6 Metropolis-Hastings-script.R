
## Example Beta(2.7, 6.3) 

{r, echo=TRUE, eval=FALSE}
## potential instumential distributions
curve(expr = dbeta(x, 2.7, 6.3), from = 0, to = 2)
curve(expr = dbeta(x, 3, 6), from = 0, to = 2, add = TRUE, col = 2, lty = 2)
curve(expr = dexp(x, rate = 3), from = 0, to = 2, add = TRUE, col = 3, lty = 2)
curve(expr = dnorm(x, mean = 0.25, 0.25), from = -2, to = 2, add = TRUE, col = 4, lty = 2)
curve(expr = dgamma(x, shape = 2.5, scale = 1/5), from = 0, to = 2, add = TRUE, col = 6, lty = 2)
legend(x = 1.2, y = 2.5, legend = c("Beta(2.7, 6.3)", "Beta(3, 6)", "Exp(1)", "N(0.25, 0.25)", "Gamma(2.5, 5)"), 
       lty = c(1, 2, 2, 2, 2), col = c(1:4, 6), merge = TRUE)




## Example Beta(2.7, 6.3) 
M = 3
par(pin = c(4.5, 3))
curve(expr = M*dgamma(x, shape = 2.5, scale = 1/5), from = 0, to = 2, ylab = "M*Gamma(2.5, 5)")
curve(expr = dbeta(x, 2.7, 6.3), from = 0, to = 2, add = TRUE, lty = 3)

#We will use $Exp(3)$ as our candidate distribution, $g$.



## Example Beta(2.7, 6.3) 
set.seed(1234)
N = 500000
## For accept-reject, we need to find a value for M

f = function(x){
  dbeta(x, 2.7, 6.3)
}

g = function(x){
  dgamma(x, shape = 2.5, scale = 1/5)
}





## Example Beta(2.7, 6.3) 
X = numeric(N)
i = 0
while(i < N){
  Y = rgamma(n = 1, shape = 2.5, scale = 1/5)
  U = runif(n = 1)
  if(U*M <= f(Y)/g(Y)){
    i = i + 1
    X[i] = Y
  }
}

qbeta(p = c(0, 0.25, 0.5, 0.75, 1), shape1 = 2.7, shape2 = 6.3) ## quantiles from Beta(2.7, 6.3)
quantile(X) ## sample mean from Accept-Reject samples




## Example Beta(2.7, 6.3) 
## see how samples from chain compare to Beta(2.7, 6.3) density
hist(X, main = "Histogram of MCMC samples", prob = TRUE, ylim = c(0, 3))
curve(expr = dbeta(x, 2.7, 6.3),  from = 0, to = 2, add = TRUE, col = "blue")


#Below is the Metropolis-Hastings implementation for this problem. 
X = numeric(N)
X[1] = rbeta(n = 1, shape1 = 2.7, shape2 = 6.3) ## initial value
for(i in 1:N){
  Y = rgamma(n = 1, shape = 2.5, scale = X[i])  #rexp(n = 1, rate = X[i])
  rho = (dbeta(x = Y, 2.7, 6.3) * dgamma(x = X[i], shape = 2.5, scale = Y) ) / 
    (dbeta(x = X[i], 2.7, 6.3) * dgamma(x = Y, shape = 2.5, scale = X[i])  )  
  
  if(runif(1) < rho){
    X[i+1] = Y
  } else{
    X[i+1] = X[i]
  }
}
qbeta(p = c(0, 0.25, 0.5, 0.75, 1), shape1 = 2.7, shape2 = 6.3) ## quantiles from Beta(2.7, 6.3)
quantile(X) ## sample mean from M-H samples


plot(X, type = "o", main = "MCMC Trace Plot", xlim = c(500,1000),
     xlab = "iterations", ylab = "X (samples obtained)")



## Example Beta(2.7, 6.3) 
## see how samples from chain compare to Beta(2.7, 6.3) density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dbeta(x, 2.7, 6.3), 
      from = 0, to = 2, add = TRUE, col = "blue")



## Example Beta(2.7, 6.3) 
#Here is a varition of the M-H algorithm used previously, except we do not let the candidate distribution depend on previous values of the chain. The candidate distribution depends only on present values of the chain, in effect $q(y | x) = q(y)$.


X = numeric(N)
X[1] = rbeta(n = 1, shape1 = 2.7, shape2 = 6.3)
for(i in 1:N){
  Y = rgamma(n = 1, shape = 2.5, scale = 1/5)
  rho = (dbeta(x = Y, 2.7, 6.3) * dgamma(x = X[i], shape = 2.5, scale = 1/5) ) / 
    (dbeta(x = X[i], 2.7, 6.3) * dgamma(x = Y, shape = 2.5, scale = 1/5) ) 
  
  if(runif(1) < rho){
    X[i+1] = Y
  } else{
    X[i+1] = X[i]
  }
  
}

quantile(X) 

## see chain transitions
plot(X, type = "o", main = "MCMC Trace Plot", xlim = c(500,1000),
     xlab = "iterations", ylab = "X (samples obtained)")

## see how samples from chain compare to Beta(2.7, 6.3) density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dbeta(x, 2.7, 6.3), from = 0, to = 2, add = TRUE, col = "blue")


#This version of the M-H algorithm is known as the **Independent Metropolis-Hastings**. This method appears a generalization of the accept-reject algorithm in the sense that the instrumental distribution is the same density $g$ as in the accept-reject algorithm. Thus, the proposed values $Y_i$ are the same, if not the accepted ones.



## Example: Gamma(4.3, 6.2) 
#Here we will compare again the Accept-Reject algorithm against the Metropolis-Hastings. Generate *N* random variables $X \sim Gamma(4.3, 6.2)$.
## For accept-reject, we need to find a value for M
## we can use `optimize` to find the maximum of our target density
maximum  = optimize(f = function(x){ dgamma(x = x, shape = 4.3, rate = 6.2)}, 
                    interval = c(0, 2), maximum = TRUE ) ## obtain maximum

M = maximum$objective
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), from = 0, 
      to = 2, col = "blue", main = "Gamma(4.3, 6.2)", xlab = "X", ylab = "Density")
abline(h = M, lty = 3, lwd = 2, col = "red")

## Example: Gamma(4.3, 6.2) 


f = function(x){
  dgamma(x = x, shape = 4.3, rate = 6.2)
}

g = function(x){
  dgamma(x = x, shape = 4, rate = 7)
}

X = numeric(N)
i = 0
while(i < N){
  Y = rgamma(n = 1, shape = 4, rate = 7)
  U = runif(1)
  if(U*M <= f(Y)/g(Y)){
    i = i + 1
    X[i] = Y
  }
}




## Example: Gamma(4.3, 6.2) 
## see how samples from chain compare to Gamma density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), 
      from = 0, to = 2, add = TRUE, col = "blue")


## Metropolis Hastings 
N = 10000
X = numeric(N)
X[1] = rgamma(n = 1, shape = 4.3, rate = 6.2)
for(i in 1:N){
  Y = rgamma(n = 1, shape = 4, rate = 7)
  rho = (dgamma(x = Y, shape = 4.3, rate = 6.2) * dgamma(x = X[i], shape = 4, rate = 7)) / 
    (dgamma(x = X[i], shape = 4.3, rate = 6.2) * dgamma(x = Y, shape = 4, rate = 7)) 
  #X[i+1] = X[i] + (Y - X[i])*(runif(1) < rho) ## equivalent to if-else statement below
  if(runif(1) < rho){
    X[i+1] = Y
  } else{
    X[i+1] = X[i]
  }
}
qgamma(p = c(0, 0.25, 0.5, 0.75, 1), shape = 4.3, rate = 6.2) #[1] 0.0000000 0.4488888 0.6405895 0.8808118       Inf
quantile(X)  ## rgamma: 0.6979356





## Example: Gamma(4.3, 6.2) 
## see chain transitions
par(mfrow = c(1,2))
plot(X, type = "o", main = "MCMC samples",
     xlim = c(1,200),
     xlab = "iterations", ylab = "X (samples obtained)")

plot(X, type = "o", main = "MCMC samples",
     xlim = c(N-200,N),
     xlab = "iterations", ylab = "X (samples obtained)")




## Example: Gamma(4.3, 6.2) 
## see how samples from chain compare to Gamma density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), 
      from = 0, to = 2, add = TRUE, col = "blue")




## Example: Gamma(4.3, 6.2) 
## potential instumential distributions
par(pin = c(4.5,3.2))
curve(expr = dgamma(x, shape = 4.3, 6.2), from = 0, to = 2, xlab = "x", ylab = "density", ylim = c(0, 1.5))
curve(expr = df(x, 4, 6), from = 0, to = 2, add = TRUE, col = 2, lty = 2)
curve(expr = dexp(x, rate = 1), from = 0, to = 2, add = TRUE, col = 3, lty = 2)
curve(expr = dgamma(x, shape = 2.5, rate = 5), from = 0, to = 2, add = TRUE, col = 6, lty = 2)
legend(x = 1.2, y = 1.5, legend = c("Gamma(4.3, 6.2)", "F(4, 6)", "Exp(1)", "N(0.25, 0.25)", "Gamma(2.5, 5)"), 
       lty = c(1, 2, 2, 2, 2), col = c(1:4, 6), merge = TRUE)


## Metropolis Hastings
X = numeric(N)
X[1] = 0.5
for(i in 1:N){
  Y = rf(n = 1, df1 = 4, df2 = 6)
  rho = (dgamma(x = Y, shape = 4.3, rate = 6.2) * df(x = X[i], df1 = 4, df2 = 6)) / 
    (dgamma(x = X[i], shape = 4.3, rate = 6.2) * df(x = Y, df1 = 4, df2 = 6)) 
  #X[i+1] = X[i] + (Y - X[i])*(runif(1) < rho) ## equivalent to if-else statement below
  if(runif(1) < rho){
    X[i+1] = Y
  } else{
    X[i+1] = X[i]
  }
}
qgamma(p = c(0, 0.25, 0.5, 0.75, 0.9), shape = 4.3, rate = 6.2)
quantile(X, probs = c(0, 0.25, 0.5, 0.75, 0.9))


## see chain transitions
plot(X, type = "o", main = "MCMC samples", xlim = c(500,1000), 
     xlab = "iterations", ylab = "X (samples obtained)")


## see how samples from chain compare to Gamma density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), from = 0, to = 2, add = TRUE, col = "blue")













# Bayesian Analysis
## O-ring Challenger Data 


failure <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 
             0, 0, 0, 1, 0, 0, 0, 0, 0)
temperature <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 
                 69, 70, 70, 70, 70, 72, 73, 75, 75,
                 76, 76, 78, 79, 81)

df = data.frame(failure, temperature)
head(df)





## O-ring Challenger Data 
fit = glm(formula = failure ~ temperature, data = df, 
          family = binomial(link = "logit"))

summary(fit)



## Output from ML estimation from Logistic Regression
## MLEs make great starting valules
a.mle <- as.numeric(fit$coefficients[1])
b.mle <- as.numeric(fit$coefficients[2])
var.a.mle <- summary(fit)$cov.scaled[1, 1]
var.b.mle <- summary(fit)$cov.scaled[2, 2]

b.hyper <- exp(a.mle + 0.577216) ## hyper parameter
## 0.577216 is "Euler’s constant"


#Let's set up the posterior and proposal distributions
# setting up functions
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



rProposal <- function(theta0){
    ## independent proposals for a and b
    #a <- log(rexp(1, 1 / b.hyper))
    a <- log(rexp(1, 1 / b.hyper)) ## log for computational purposes
    b <- rnorm(1, b.mle, sqrt(var.b.mle))
    return(c(a, b))
}




## O-ring Challenger Data 

Now we run the M-H algorithm


## Metropolis-Hastings set up to run
# Run Metropolis-Hastings
N = 1000000
BurnIn =  5000

x0 <- c(a.mle, b.mle) ## x0: initializing values for MH algorithm
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


## O-ring Challenger Data | Summary 



summary(MH.Results) ## summary of parameters


par(mfrow = c(1,2))
hist(alphaMH, freq=FALSE, col="gray", border="white", xlab=expression(alpha), main = "")
plot(alphaMH, type="l", col="gray", xlab="Iteration", ylab=expression(alpha), main = "")
lines(1:N, cumsum(alphaMH) / (1:N))



## Exercise Student's t density with v degrees of freedom
#Calculate the mean of a t distribution with $v = 4$ degrees of freedom using a M-H algorithm with candidate densities $N(0,1)$ and $t_{v = 2}$.

# Appendix

## Solution to Student's t density with v degrees of freedom

#Calculate the mean of a t distribution with $v = 4$ degrees of freedom using a M-H algorithm with candidate densities $N(0,1)$ and $t_{v = 2}$.



set.seed(987)
N = 10^6

#dt(x = x, df = 4)

X = numeric(N)
X[1] = rnorm(1) ## initialize the starting value

for(i in 1:N){
    Y = rnorm(1) ## independent of X_i
    rho = (dt(Y, df = 4) * dnorm(X[i])) /
            (dt(X[i], df = 4) * dnorm(Y))
    U = runif(1)
    if(U <= rho){
        X[i+1] = Y
    } else{
        X[i+1] = X[i]
    }

}


plot(density(X), type = "l", 
        lty = 2, main = "M-H with N(0,1) candidate")
curve(dt(x, df = 4), add = TRUE, col = 4)




