## Monte Carlo Integration 

# If we now simulate this, we will see we approximate the true solution $\int_0^2h(x) dx = 8$
n = 50000 ## sample size
h <- function(x) { 3*x^2 } ## function of interest, h(x)
X <- runif(n = n, min = 0, max = 2)  ## samples from f(x)
v = 2 * mean(h(X)) ## 2 * E[h(x)]
print(v) ## approximately 8
integrate(f = h, lower = 0, upper = 2)




N = 10000 ## sample size
h <- function(x) { 3*x^2 } ## function of interest, h(x)
X <- runif(n = N, min = 0, max = 2)  ## samples from f(x)
h_values <- 2 * h(X)

cumMean <- function(x, n){
  num = cumsum(x)  ## numerator
  denom = 1:n  ## denominator
  result = num/denom
  return(result)
}







cumSE <- function(x, n){
  m = mean(x)
  num = sqrt(cumsum((x - m)^2)) ## numerator
  denom = 1:n ## denominator
  result = num/denom ## cummulative mean of (x_i - theta)^2
  return(result)
}

thetas = cumMean(h_values, N)
SE = cumSE(h_values, N)







plot(x = 1:N, y = thetas, type = "l", ylim = c(4, 12), xlab = "Number of Samples",
     ylab = "Estimate of theta",
     main = "Estimate of Theta with 95% CI")
lines(x = 1:N, y = thetas + 1.96*SE, col = "gray") ## CI
lines(x = 1:N, y = thetas - 1.96*SE, col = "gray") ## CI






## final estimate
thetaHat = mean(h_values)
se <- sd(x = h_values)/sqrt(N)
ci <- thetaHat + 1.96*c(-1,1) * se

print(thetaHat) ## theta estimate
print(ci) ## 95% CI




x <- seq(from = -10, to = 10, by = 0.1)

out = dnorm(x)
out[x > -1] = 0

par(mfrow = c(1,2))
plot(x, dnorm(x), type = "l", main = "Phi(X)", ylab = "density")
polygon(x, out, col = "blue")


out = dnorm(x)
out[x < 1] = 0

plot(x, dnorm(x), type = "l", main = "Phi(X) = 1 - Phi(-X)", , ylab = "density")
polygon(x, out, col = "red")




N = 10000
x = 2.5
t = runif(N, 0, 2.5)
p = 0.5 + 1/sqrt(2*pi) * mean(x * exp(-t^2 / 2))

print(p)
print(x = pnorm(q = 2.5, lower.tail = TRUE))





N = 10000
x = 3.2
Z = rnorm(N)
p_val = mean(Z >= x)
print(p_val)
print(sum(Z > x)) ## few samples with Z > x






##  Calculating tail probabilities 
set.seed(9124555) ## won't show up

N = 10000; x = 4
t = rnorm(N) >= x ## looking at the lower tail
p_val = sum(t)/N
print(sum(t)) ## only 1 sample greater or equal to 4
print(p_val)

## compare this with
RPval = pnorm(q = 4, lower.tail = FALSE) ## lower tail
print(RPval)







N = 10^4 
t = 4
X = rexp(N) + t
w = dnorm(X)/dexp(X - t) ## importance sampling weights dnorm/dexp(x-4)

theta = mean(w)
cumTheta = cumsum(w) / 1:N


plot(cumTheta, type = "l", main = "N(0,1) tail probabilities",
     xlab = "Number of samples", ylab = "p(X>=4)",
     ylim = c(2.5 * 10^(-5), 4 * 10^(-5)))
abline(a = pnorm(t, lower.tail = FALSE), b = 0, col = "blue")

## Exercise: Compare the results
#print(pnorm(t, lower.tail = FALSE))
#print(theta) 


plot(density(X, weights=w/sum(w)), main = "Density of samples * weights",
     xlab = "X", ylab = "Density",
     xlim = c(0,16)) ## plot density of samples from Exp(x - 4) * weights







##  Calculating Expectation of Laplace Distribution  
f <- function(x){
  out = (1/2) * exp(x = -1 * abs(x))
  return(out)
}

curve(expr = f, from = -10, to = 10, type = "l",
      main = "Laplace Distribution")

curve(expr = dnorm(x, sd = 1), from = -10, to = 10, type = "l",
      lty = 2, col = 2,
      add = TRUE)

curve(expr = dnorm(x, sd = 2), from = -10, to = 10, type = "l",
      lty = 3, col = 3,
      add = TRUE)


curve(expr = dnorm(x, sd = 4), from = -10, to = 10, type = "l",
      lty = 4, col = 4,
      add = TRUE)

legend(7, 0.4, c("sd = 1", "sd = 2", "sd = 4"), col = c(2, 3, 64),
       text.col = "gray30", lty = c(2, 3, 4), pch = c(NA, NA, NA),
       merge = TRUE)


N = 10000
x = rnorm(N)
w = f(x)/dnorm(x, sd = 1) 
theta <- mean(x = x^2 * w)
print(theta)


x = rnorm(N, sd = 2)
w = f(x)/dnorm(x, sd = 2) 
theta <- mean(x = x^2 * w)
print(theta)




x = rnorm(N, sd = 4)
w = f(x)/dnorm(x, sd = 4)
theta <- mean(x = x^2 * w)
print(theta)








N = 10^4
f <- function(x){
  out = (1/2) * exp(x = -1 * abs(x))
  return(out)
}
pos_neg <- sample(x = c(-1,1), size = N, replace = TRUE)
x = rexp(n = N)
double_exp = pos_neg * x 

hist(double_exp, probability = TRUE, ylim = c(0,0.7),
     main = "Double Exponential from Exp(1)",
     xlab = "Double Exponential", ylab = "Density")
lines(density(double_exp))
curve(expr = f, from = -4, to = 4, add = TRUE, col = 4, lty = 4)


theta = mean(double_exp^2)
print(paste("theta estimate: ", theta))






## Tail Probabilities of a *t*-distribution 

n = 5*10^4
x = rexp(n = n, rate = 1)
alpha = 4
y = x + alpha
par(mfrow = c(1,2), mar = c(4,4,2,2))

hist(x, probability = TRUE, xlim = c(0, 14), col = "gray", border = "white")
hist(y, probability = TRUE, xlim = c(0, 14), col = "gray", border = "white")

pvals = cumsum(w/1:n)
sterr = sqrt(cumsum((w - pvals)^2))/(1:n)

par(mar = c(4,4,2,2))
plot(x = 1:n, y = pvals, type = "l", lty = 2, 
     xlab = "samples", ylab = "p-val estimates", ylim = c(0, 0.001))
lines(pvals + 1.96*sterr, col = "blue", lty = 3)
lines(pvals - 1.96*sterr, col ="blue", lty = 3)











## Normal-Cauchy Bayes estimator 
N = 10^5
X = 2.5

thetas = rnorm(N, X, sd = 1)
weights = dcauchy(x = thetas, location = 0, scale = 1) ## 1/(pi * (1 + s^2))
sumOfWeights = sum(weights)

posteriorDist =  (weights )/(sumOfWeights)
hist(posteriorDist, probability = TRUE, 
     col = "gray", border = "white",
    main = "Posterior Distribution of thetas",
    xlab = "theta", ylab = "density")



## The expectation of the posterior distribution
posteriorThetaEst = sum((thetas * weights)/ sumOfWeights)
print(posteriorThetaEst)


## Effective Sample Size
ESS = 1 / sum( (weights / sumOfWeights)^2) ## Effective Sample Size
print(ESS)




posteriorTheta = cumsum(thetas * weights) / cumsum(weights)

plot(posteriorTheta, type = "l", main = "posterior theta estimates",
     ylim = c(1.6, 1.85), col = "gray")



X = 2.5
N = 10^4
M = 1000
thetas = matrix(data = rnorm(N*M), ncol = M) + X ## N(0,1) + 2.5
weights = dcauchy(x = thetas, location = 0, scale = 1)
posteriorTheta = apply(thetas * weights, 2, cumsum) / apply(weights, 2, cumsum)



## look at the first column estimates
plot(posteriorTheta[,1], main = "posterior theta estimates",
    ylim = c(1, 2.5), type = "l")
## grab the 2.5 and 97.5 quantiles for a 95% CI
confInts = apply(posteriorTheta, 1, quantile, c(0.025, 0.975)) 
polygon(x = c(1:N, N:1), y = c(confInts[1,], rev(confInts[2,])), border = "gray")




