## importance sampling example
## Normal with Cauchy prior


N = 10^5
X = 2.5

thetas = rnorm(N, X, sd = 1)

weights = dcauchy(x = thetas, location = 0, scale = 1) ## 1/(pi * (1 + s^2))
sumOfWeights = sum(weights)

posteriorDist =  (weights )/(sumOfWeights)
hist(posteriorDist, probability = TRUE, 
    main = "Posterior Distribution of thetas",
    xlab = "theta", ylab = "density")


## The expectation of the posterior distribution
posteriorThetaEst = sum((thetas * weights)/ sumOfWeights)
print(posteriorThetaEst)


## Effective Sample Size
ESS = 1 / sum( (weights / sumOfWeights)^2)
print("Effective Sample Size")
print(ESS)


posteriorTheta = cumsum(thetas * weights) / cumsum(weights)

plot(posteriorTheta, type = "l", main = "posterior theta estimates",
    ylim = c(1.6, 1.85))


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
polygon(x = c(1:N, N:1), y = c(confInts[1,], rev(confInts[2,])))




