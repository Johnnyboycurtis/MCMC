X <- runif(100000,0,10) 
Y <- 10*exp(-2*abs(X-5))
c( mean(Y), var(Y) ) 
w <- function(x) dunif(x, 0, 10)/dnorm(x, mean=5, sd=1)
f <- function(x) 10*exp(-2*abs(x-5))
X=rnorm(1e5,mean=5,sd=1)
Y=w(X)*f(X)
c( mean(Y), var(Y) ) 
X <- rnorm(1e5, sd=2) 
Y <- (X^2) * .5 * exp(-abs(X))/dnorm(X, sd=2) 
mean(Y)
IS <- function(X, C, res) 
{
 n <- length(X) 
 log.posterior <- function(t) (sum(X)+4)*log(t) + (n*(10-mean(X))-2)*log(1-t)
 log.g <- function(t) dbeta(t,C*mean(X),C*(10-mean(X)), log=TRUE)
 log.w <- function(t) log.posterior(t) - log.g(t)
 U <- rbeta(res, C*mean(X), C*(10-mean(X)))
 LP <-  log.w(U) 
 
    # factor out the largest value to prevent numerical underflow
    w <- max(LP)
    LP <- LP - w 
 I <- mean( exp(LP)*U )/mean(exp(LP))
 ESS <- mean( ( (exp(LP)/mean(exp(LP))) - 1)^2 )
 Z = exp(LP)/mean(exp(LP))
 sig.sq <- (1/res)*sum( (Z-I)^2 ) 
 se <- sqrt( sig.sq/res ) 
 return(c(I - 1.96*se, I + 1.96*se, ESS))
}
X = rbinom(100, 10, .4)
const <- seq(.1, 200, length=2000) 
A <- matrix(0, 2000, 2) 
for(j in 1:2000)
{
 
 OUT <- IS(X, const[j], 1000) 
 ESS <- OUT[3]
 V <- (OUT[2]-OUT[1])/3.92
 A[j,] <- c(ESS, V)
} 
plot(const, A[,1], ylab="ESS", xlab="c", main="ESS vs. c", col=2, type="l")
abline(h=5)
plot(const, A[,2], xlab="c", ylab="Standard error estimate", main="se vs. c", col=4, type="l")
IS(X, const[which.min(A[,2])], 1000)[1:2]
IS <- function(X, C, res) 
{
 n <- length(X) 
 log.posterior <- function(t) ( (4-n)/2 )*log(t) - (1/(2*t)) * (t^2 + sum(X^2))
 a = C*var(X); b = C; 
 log.g <- function(t) dgamma(t,a,b,log=TRUE) 
 log.w <- function(t) log.posterior(t) - log.g(t)
 U <- rgamma(res, a, b)
 LP <-  log.w(U) 
    w <- max(LP)
    LP <- LP - w 
 I <- mean( exp(LP)*U )/mean(exp(LP))
 ESS <- mean( ( (exp(LP)/mean(exp(LP))) - 1)^2 )
 Z = exp(LP)/mean(exp(LP))
 sig.sq <- (1/res)*sum( (Z-I)^2 ) 
 se <- sqrt( sig.sq/res ) 
 return(c(I - 1.96*se, I + 1.96*se, ESS))
}
X = rnorm(100,mean=0,sd=2)
const <- seq(.05, 20, length=500) 
A <- matrix(0, 500, 2) 
for(j in 1:500)
{
 
 OUT <- IS(X, const[j], 1000) 
 ESS <- OUT[3]
 V <- (OUT[2]-OUT[1])/3.92
 A[j,] <- c(ESS, V)
} 
plot(const, A[,1], ylab="ESS", xlab="c", main="ESS vs. c", col=2, type="l")
abline(h=5)
plot(const, A[,2], xlab="c", ylab="standard error vs. c", col=4, type="l")
IS(X, const[which.min(A[,2])], 1000)[1:2]


