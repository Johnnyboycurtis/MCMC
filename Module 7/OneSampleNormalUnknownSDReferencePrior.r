library(coda)
source("Bayesfunctions.r") # You must use the correct working directory

n <-  50
tau <- 1 /13^2
outcomes <- floor(rnorm(n, 110, 1/sqrt(tau)))
# the value of mu used to simulate the data is 110 and the value of sigma is 13.

summary(outcomes)
hist(outcomes)

y.bar <- mean(outcomes)  # MLE

# Assume now that the standard deviation (precision) of the y's is unknown

# We need a prior on theta=(mu, tau)

# We will use an "improper" reference prior on theta=(mu, tau)

# f(mu, tau) = 1 / tau ; this is a density function but not a probability density function since it doesn't integrate to 1.

# We sample from the posterior of (mu, tau) using sequential sampling

s <- sd(outcomes)

n.iterations <-1000000
tau.posterior.sample <- rep(NA, n.iterations)
mu.posterior.sample <- rep(NA, n.iterations)

for(i in 1:n.iterations){
tau.posterior.sample[i] <- rgamma(1, shape = (n-1) / 2, rate = (n-1) * s^2 /2)
mu.posterior.sample[i] <- rnorm(1, mean = y.bar, sd = 1 / sqrt(n * tau.posterior.sample[i]))
}

summary.Bayes(mu.posterior.sample, HPD.p=0.95)
summary.Bayes(tau.posterior.sample, HPD.p=0.95)

# We care about the standard deviation sigma
sigma.posterior.sample <- 1 / sqrt(tau.posterior.sample)
summary.Bayes(sigma.posterior.sample, HPD.p=0.95)


# Exact 95% PI for mu
y.bar + qt(0.975, df=n-1) * s / sqrt(n) * c(-1, 1)
# Estimated 95% PI for mu
HPDinterval(as.mcmc(mu.posterior.sample))


par(mfrow=c(1, 2))
#marginal posterior for mu
plot(density(mu.posterior.sample), type="l", main=expression(paste("Posterior Density of ", mu)), xlab=expression(mu), ylab="density")
lines(HPDinterval(as.mcmc(mu.posterior.sample)), c(0,0), col="red", cex=2, lty=3)  # add PI line
legend("topright", legend=c("posterior", "95% HPD interval"), lty=c(1, 2), col=c("black", "red"))


# Exact 95% PI for sigma (not necessarily the HPD)
lower <- sqrt((n-1) / qchisq(0.975, df=n-1)) * s
upper <- sqrt((n-1) / qchisq(0.025, df=n-1)) * s
c(lower, upper)

# Estimated 95% HPD interval for sigma
HPDinterval(as.mcmc(sigma.posterior.sample))

# marginal posterior for sigma
plot(density(sigma.posterior.sample), type="l", main=expression(paste("Posterior Density of ", sigma)), xlab=expression(sigma), ylab="density")
lines(HPDinterval(as.mcmc(sigma.posterior.sample)), c(0,0), col="red", cex=2, lty=3)  # add PI line
legend("topright", legend=c("posterior", "95% HPD interval"), lty=c(1, 2), col=c("black", "red"))


# Exact 95% PI for the new observation
y.bar + qt(0.975, df=n-1) * s * sqrt(1 + 1/n) * c(-1, 1)

# Estimated 95% PI for the new observation
HPDinterval(as.mcmc(y.tilde.sample))


# Predictive distribution
y.tilde.sample <- y.bar + s * sqrt(1 + 1/n) * rt(n.iterations, df=n-1) # random sample from the predictive density
par(mfrow=c(1, 1))
plot(density(y.tilde.sample), main="Predictive Density", xlab="new observation", ylab="density")
lines(HPDinterval(as.mcmc(y.tilde.sample)), c(0,0), col="red", cex=2, lty=3)  # add Prediction Interval line
legend("topright", legend=c("predictive density", "95% Prediction Interval"), lty=c(1, 2), col=c("black", "red"))

# Summary of the predictive distribution
summary(y.tilde.sample)





