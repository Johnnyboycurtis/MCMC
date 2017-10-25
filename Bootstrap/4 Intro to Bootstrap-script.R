
#We will use the `mtcars` data set to illustrate a simple implementation. 
data("mtcars")
mpg = mtcars$mpg
n = length(mpg)
print(mean(mpg))
hist(x = mpg, probability = TRUE, xlab = "MPG", main = "Histogram of MPG")



## Bootstrap Example

B = 1000 ## number of bootstraps
results = numeric(B) ## vector to hold results
for(b in 1:B){
  i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
  bootSample = mpg[i] ## get data
  thetaHat = mean(bootSample) ## calculate the mean for bootstrap sample
  results[b] = thetaHat ## store results
}

hist(x = results, probability = TRUE, 
     main = "Bootstrapped Samples of Mean_mpg",
     xlab = "theta estimates")





library(ggplot2, quietly = TRUE) ## for graphics
mtcars$am <- as.factor(mtcars$am) ## Transmission (0 = automatic, 1 = manual
fit = lm(formula = mpg ~ wt + am, data = mtcars)
data.frame(coefficients = coefficients(fit), CI = confint(fit), check.names = FALSE)



## Paired Bootstrapping 


qplot(x = as.factor(am), y = mpg, data = mtcars, geom = "boxplot",
      main = "Boxplot: MPG ~ AM", ylab = "MPG", xlab = "AM",
      colour = am)



## Paired Bootstrapping 


qplot(x = wt, y = mpg, data = mtcars, geom = c("point", "smooth"),
      main = "Boxplot: MPG ~ Weight", ylab = "MPG", xlab = "Weight",
      method = "lm", formula = y~x)


## save coefficients
beta_int = coefficients(fit)[1]
beta_wt = coefficients(fit)[2]
beta_am = coefficients(fit)[3]

n = dim(mtcars)[1] ## number of obs in data
B = 1000 ## number of bootstrap samples

results = matrix(data = NA, nrow = B, ncol = 3, 
                 dimnames = list(NULL, c("Intercept", "wt", "am")))

## begin bootstrap for-loop
for(b in 1:B){
  i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
  temp = mtcars[i,] ## temp data set
  temp_model =  lm(formula = mpg ~ wt + am, data = temp) ## train model
  coeff = matrix(data = coefficients(temp_model), ncol = 3) ## get coefficients
  results[b,] = coeff ## save coefficients in matrix
}





## Paired Bootstrapping 




results <- data.frame(results, check.names = FALSE)

summary(results) ## take a look at the samples

boot_int = results[,"Intercept"]
boot_wt = results[,"wt"]
boot_am = results[,"am"]




## Paired Bootstrapping 
par(mfrow = c(2,2))
hist(boot_int, main = "Bootstrapped Coefficients for Intercept",
     xlab = "Coefficients for Intercept", probability = TRUE)
abline(v = coefficients(fit)[1], col = "black", lty=2)


hist(boot_wt, main = "Bootstrapped Coefficients for Weight",
     xlab = "Coefficients for Weight", probability = TRUE)
abline(v = coefficients(fit)[2], col = "blue", lty=2)


hist(boot_am, main = "Bootstrapped Coefficients for AM = 1",
     xlab = "Coefficients for Automatic Transmission", probability = TRUE)
abline(v = coefficients(fit)[3], col = "green", lty=2)





## Paired Bootstrapping 
#Now we can estimate bias for each parameter estimate. Define Bias as $Bias(\theta) = E[\theta^*] - \theta$, where in our scenario we have $Bias(\hat\theta) = E[\hat\theta^*] - \hat\theta$. Our bootstrap bias corrected estimates are then $\hat\theta_{BC} = \hat\theta - Bias(\hat\theta)$.

bias_int = mean(boot_int - beta_int)
print(bias_int)

bias_wt = mean(boot_wt - beta_wt)
print(bias_wt)

bias_am = mean(boot_am - beta_am)
print(bias_am)


## incorportate our bias into the coefficients
## we now have bias corrected coefficients

intercept = beta_int - bias_int
print(intercept)

wt = beta_wt - bias_wt
print(wt)

am = beta_am - bias_am
print(am)


