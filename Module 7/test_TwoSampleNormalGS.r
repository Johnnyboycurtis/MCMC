## Two-Sample Normal Model with Proper Independent Priors : R-Code
library(coda)
#setwd("C:/Users/Jonathan/OneDrive/SDSU/Stat 676 Bayesian Statistics/R Code")
source("Bayesfunctions.r")


## Diasorin as a diagnostic test for low bone turnover in kidney dialysis patients
## Comparing two means
## Ch 5, Section 2

low <- c(91, 46,  95,  60, 33, 410, 105, 43, 189,1097, 54,178, 114, 137, 233, 101, 25,70,357)
normal <- c(370, 267, 99,157, 75,1281, 48, 298, 268, 62,804,430,171,694,404)
n <- c(length(low), length(normal))
diasorin <- c(low, normal)
group <- c(rep(1, n[1]), rep(2, n[2]))
group <- factor(group, 1:2, c('Low','Normal'))

####################################
#     Normality. Fig 5.2           #
####################################

par(mfrow=c(2,2), pty="s")
boxplot(diasorin~group,range=0, ylab="Diasorin", lwd=1.5)
boxplot(log(diasorin)~group,range=0, ylab="log(Diasorin)", lwd=1.5)

qqnorm(log(diasorin[group=="Low"]),main = "Low Group \n Log transformed data", lwd=2)
qqline(log(diasorin[group=="Low"]), lwd=2)

qqnorm(log(diasorin[group=="Normal"]),main = "Normal Group \n Log transformed data", lwd=2)
qqline(log(diasorin[group=="Normal"]), lwd=2)


dev.off() # To go back to the plot default settings

####################################
#     Descriptive Statistics       #
####################################

# install.packages("Hmisc")
library(Hmisc)

g.summary <- function(x){
  return(c(MEAN=mean(x,na.rm=TRUE),SD=sd(x,na.rm=TRUE),MEDIAN=median(x,na.rm=TRUE),MIN=min(x,na.rm=TRUE),MAX=max(x,na.rm=TRUE),N=length(x)))
}

summarize(diasorin,by=group, FUN=g.summary)
summarize(log(diasorin),by=group, FUN=g.summary)


# Summary stats to be used in the Gibbs sampler
y.bar <- c(mean(log(low)), mean(log(normal)))
s <- c(sd(log(low)), sd(log(normal)))


####################################
#     Prior Densities      #
####################################

# Hyper-parameters
# mu=(mu1, mu2) prior
find.normal(prior.mean=log(130), percentile=log(142), p=0.95)
find.normal(prior.mean=log(220), percentile=log(240), p=0.95)

mu.prior.mean <- c(4.87,  5.39)
mu.prior.sd <- c(0.0537, 0.0529)
mu.prior.precision <- 1 / mu.prior.sd^2

par(mfrow=c(1,2))
for(i in 1:2){
plot(density(rnorm(100000, mean=mu.prior.mean[i], sd=mu.prior.sd[i])),
     substitute(paste("Prior Density of ", mu[a]), list(a=levels(group)[i])),
     xlab=substitute(mu[a], list(a=levels(group)[i])),
     ylab="density", xlim=c(4.6, 5.6))
}

# tau = (tau1, tau2) and sigma = (sigma1, sigma2) prior
normal.percentile.to.sd(mean.value=log(130), percentile=log(170), p=0.95) # Returns prior mode guess for sigma_1
normal.percentile.to.sd(mean.value=log(130), percentile=log(200), p=0.95) # Returns prior percentile guess for sigma_1

gamma.parameters1 <- find.tau.gamma(prior.sigma.mode=0.163, sigma.percentile=0.262, p=0.95) # Returns shape and rate parameters for the Gamma distribution of tau
gamma.parameters1

tau.prior.shape <- c(NA, NA)
tau.prior.rate <- c(NA, NA)

tau.prior.shape[1] <- gamma.parameters1$a
tau.prior.rate[1] <-  gamma.parameters1$b

normal.percentile.to.sd(mean.value=log(220), percentile=log(280), p=0.95) # Returns prior mode guess for sigma_2
normal.percentile.to.sd(mean.value=log(220), percentile=log(300), p=0.95) # Returns prior percentile guess for sigma_2

gamma.parameters2 <- find.tau.gamma(prior.sigma.mode=0.1466161, sigma.percentile=0.1885608, p=0.95) # Returns shape and rate parameters for the Gamma distribution of tau
gamma.parameters2

tau.prior.shape[2] <- gamma.parameters2$a
tau.prior.rate[2] <-  gamma.parameters2$b

tau.prior.shape
tau.prior.rate

par(mfrow=c(2, 2))
for(i in 1:2){
plot(density(rgamma(100000, shape=tau.prior.shape[i], rate=tau.prior.rate[i])), main=substitute(paste("Prior Density of ", tau[a]), list(a=levels(group)[i])), xlab=substitute(tau[a], list(a=levels(group)[i])), ylab="density", xlim=c(0, 100))
plot(density(1/sqrt(rgamma(100000, shape=tau.prior.shape[i], rate=tau.prior.rate[i]))), main=substitute(paste("Prior Density of ", sigma[a]), list(a=levels(group)[i])), xlab=substitute(sigma[a], list(a=levels(group)[i])), ylab="density", xlim=c(0.07, 0.4))
}

####################################
#    Gibbs Sampler     #
####################################
burn.in <- 10000
net.iterations <- 100000
n.iterations <- burn.in + net.iterations
tau.GS <- matrix(NA, nrow=n.iterations, ncol=2) # We now store the values in a n.iterations x 2 matrix
mu.GS <- matrix(NA, nrow=n.iterations, ncol=2)

# Starting values

tau.GS[1, ] <- c(1, 1)
mu.GS[1, ] <- y.bar # this is a vector!

tau.post.shape <- tau.prior.shape+ 0.5 * n

start.time <- Sys.time()
for(i in 1:(n.iterations - 1)){
  for(j in 1:2){
  tau.GS[i+1, j] <- rgamma(1, shape = tau.post.shape[j] , rate = tau.prior.rate[j] + 0.5 * (n[j] * (y.bar[j] - mu.GS[i, j])^2 + (n[j]-1) * s[j]^2))
  weight <- n[j] * tau.GS[i+1, j] / (n[j] * tau.GS[i+1, j] + mu.prior.precision[j])
  mu.GS[i+1, j] <- rnorm(1, mean = weight * y.bar[j] + (1-weight) * mu.prior.mean[j], sd = 1 / sqrt(n[j] * tau.GS[i+1, j] + mu.prior.precision[j]))
  }
}
Sys.time() - start.time # how long did it take to run the Gibbs sampler?

# We care about the standard deviation sigma
sigma.GS <- 1 / sqrt(tau.GS)

mu.GS <- mu.GS[(burn.in+1):n.iterations, ]
tau.GS <- tau.GS[(burn.in+1):n.iterations, ]
sigma.GS <- sigma.GS[(burn.in+1):n.iterations, ]

apply(mu.GS, MARGIN=2, FUN=summary.Bayes) # apply the summary.Bayes function to each column (group)
apply(sigma.GS, MARGIN=2, FUN=summary.Bayes)

# Summarize the mean difference
mu.diff.GS <- mu.GS[, 2] - mu.GS[, 1]
summary.Bayes(mu.diff.GS)

# Summarize the sd ratio
sigma.ratio.GS <- sigma.GS[, 2] / sigma.GS[, 1]
summary.Bayes(sigma.ratio.GS)


####################################
#   Estimated Posterior Densities  #
####################################

dev.off()  # return to default plot settings

par(mfrow=c(2,1))
#marginal posterior for each mu_j
for(j in 1:2){
plot(density(mu.GS[, j]), type="l", main=substitute(paste("Posterior Density of ", mu[a]), list(a=levels(group)[j])), xlab=substitute(mu[a], list(a=levels(group)[j])), ylab="density", xlim=range(mu.GS))
lines(HPDinterval(as.mcmc(mu.GS[,j])), c(0,0), col="red", cex=2, lty=3)  # add PI line
legend("top", legend="95% HPD interval", lty=2, col="red", cex=0.4)
}

par(mfrow=c(2,1))
# marginal posterior for each sigma_j
for(j in 1:2){
plot(density(sigma.GS[, j]), type="l", main=substitute(paste("Posterior Density of ", sigma[a]), list(a=levels(group)[j])), xlab=substitute(sigma[a], list(a=levels(group)[j])), ylab="density", xlim=range(sigma.GS))
lines(HPDinterval(as.mcmc(sigma.GS[,j])), c(0,0), col="red", cex=2, lty=3)  # add PI line
legend("topright", legend="95% HPD interval", lty=2, col="red", cex=0.4)
}



