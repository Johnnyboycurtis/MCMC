getwd()  # Ask for current working directory
setwd("~/Dropbox/SDSU/SDSU_courses/STAT676/R-code") # Set working directory ; Use your own path
getwd()  # Verify working directory

library(nlme)
library(coda)
library(rjags)



# Read the cleanPIMA data set
mydata<-read.table("./Binomial Models/PIMA.csv", header=TRUE)
attach(mydata)
# Summary of data
summary(mydata)

# Check level order
levels(test)

# Logistic Regression
fit.glucose <- glm( test ~ glucose, family=binomial("logit"))

# Summary of fit
summary(fit.glucose)

# Obtain predicted probabilities
predicted<-predict(fit.glucose, type="response")

# Plot observations predicted probabilities
plot(glucose[complete.cases(glucose)], as.numeric(test)[complete.cases(glucose)] - 1,
     xlab="Glucose level", ylab="Probability of testing positive for diabetes")
# Add predicted probabilities to the plot
points(glucose[complete.cases(glucose)], predicted)

#Estimate of probabilities for new observations
newdata <- data.frame(glucose=c(100, 125, 150, 175, 200))
predict(fit.glucose, newdata, type="response")


# Confidence intervals for regression coefficients
confint(fit.glucose)
# Confidence intervals for odds ratio (comparing groups of people whose glucose value differs by 1 unit)
exp( cbind( OR=coef(fit.glucose), confint(fit.glucose) ) )



#######################################
#      MCMC Using JAGS                #
#######################################

# We have to take out the people missing the TEST or GLUCOSE variables

test01 <- as.numeric(test) - 1  # Jags wants 0s and 1s and not negatives and positives
complete.data <- na.omit((data.frame(test01, glucose)))

summary(complete.data)

# Create file with model specification
cat( "
     var
     test[N], glucose[N], beta[2], OR ;
     model{
     # Likelihood specification
     for (i in 1:N)
     {
     test[i] ~ dbern(theta[i])
     logit(theta[i]) <- beta[1] + beta[2] * glucose[i]
     }
     # Prior specification
     for(i in 1:2){
     beta[i] ~ dnorm(0.0, 1E-6)  # Reference prior
     }
     OR <- exp(beta[2]) # parameter of interest
     }",
file="LogisticRegressionPimaExampleReferencePrior.jags")


jagsfit <- jags.model(file = "LogisticRegressionPimaExampleReferencePrior.jags",
                      data = list('test' = complete.data$test,
                                  'glucose' = complete.data$glucose,
                                  'N' = length(complete.data$test)),
                       inits = list('beta'= c(0, 0)),
                      n.chains=3,
                      n.adapt=0
)

start.time <- Sys.time()
update(jagsfit, 10000)
Sys.time()  - start.time


start.time <- Sys.time()
MCMC.out <- coda.samples(jagsfit,
                         var = c("beta", "OR"),
                         n.iter = 10000,
                         thin = 1)
Sys.time()  - start.time
summary(MCMC.out)


######## Analysis with informative priors

X.tilde <- matrix(c(1, 1, 100, 200), nrow=2, ncol=2)
X.tilde.inv <- solve(X.tilde)

a <- c(5 , 8) # prior successes
b <- c(20 , 4) # prior failures

# Create file with model specification
cat( "
     var
     test[N], glucose[N], beta[2], OR, tildetheta[2], X.tilde.inv[2, 2], a[2], b[2] ;
     model{
     # Likelihood specification
     for (i in 1:N)
     {
     test[i] ~ dbern(theta[i])
     logit(theta[i]) <- beta[1] + beta[2] * glucose[i]
     }
     # Prior specification
     for(i in 1:2){
     tildetheta[i] ~ dbeta(a[i], b[i])  # Informative prior
     }
     beta[] <- X.tilde.inv[ , ] %*% logit(tildetheta[])
     OR <- exp(beta[2]) # parameter of interest
     }",
     file="LogisticRegressionPimaExampleInformativePrior.jags")


jagsfit <- jags.model(file = "LogisticRegressionPimaExampleInformativePrior.jags",
                      data = list('test' = complete.data$test,
                                  'glucose' = complete.data$glucose,
                                  'N' = length(complete.data$test),
                                  "a" = a,
                                  'b' = b,
                                  'X.tilde.inv' = X.tilde.inv),
                      inits = list('tildetheta'= c(0.5, 0.5)),
                      n.chains=3,
                      n.adapt=0
)

start.time <- Sys.time()
update(jagsfit, 100000)
Sys.time()  - start.time


start.time <- Sys.time()
MCMC.I <- coda.samples(jagsfit,
                         var = c("beta", "OR"),
                         n.iter = 100000,
                         thin = 1)
Sys.time()  - start.time
summary(MCMC.I)




###############################################################
#                 MCMC Diagnostics                            #
###############################################################

# library(coda)

gelman.diag(MCMC.I, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)


png('GewekePlotLogistic.png') # Save plot as .png file
geweke.plot(MCMC.I[[1]], frac1 = 0.1, frac2 = 0.5, nbins = 20,
            pvalue = 0.05, auto.layout = TRUE)
dev.off()


OR <- MCMC.I[ , 1, drop=FALSE] # select chains for first variable
png('ORMCMCI.png')
plot(OR)
dev.off()

# As an exercise, look at the MCMC diagnostics for the model with the non-informative prior.
