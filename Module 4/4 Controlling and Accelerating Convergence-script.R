## Monitoring Convergence 
#As a toy example, consider the simple function $h(x) = [cos(50x) + sine(50x)]^2$. Using a simple MC algorithm, we can estimate $\theta = \int_0^1 h(x)$. Let us generate $n$ samples $x_1, ..., x_n \sim Unif(0,1)$, such that $\displaystyle \theta = E[h(x)] \approx \frac{1}{n} \sum^n h(x_i)$. 


set.seed(3456)
n = 80000
x = runif(n)
h = function(x){
  v = (cos(50*x) + sin(50*x))^2
  return(v)
}

#thetaHat = mean(h(x)) # theta hat

theta_est = cumsum(h(x))/1:n ## cumulative mean

se = sqrt( cumsum((h(x) - theta_est)^2 ) / 1:n ) / sqrt(1:n)



## Monitoring Convergence 
par(pin = c(5,3))
plot(x = 1:n, y = theta_est, type = "l", lty=2, ylim = c(0.7, 1.2), xlim = c(0, 10000))
lines(x = 1:n, y = theta_est - 1.96*se, col = "blue")
lines(x = 1:n, y = theta_est + 1.96*se, col = "blue")
legend(x = 6000, y = 0.85, legend = c(expression(hat(theta)), "CI"), 
       border = "white", col = c("black", "blue"), lty = c(2, 1))






## While running, monitor your computer's resource manager


## parallel monte carlo samples
M = 200L
X = matrix(data = runif(n*M), nrow = n, ncol = M)
h_samples = h(X)
thetaEstimates = apply(X = h_samples, MARGIN = 2, FUN = function(v){ cumsum(v)/1:n } )

parallelCI = t(apply(X = thetaEstimates, MARGIN = 1, FUN = quantile, c(0.025, 0.50, 0.975)))

summary(parallelCI)

integrate(f = h, lower = 0, upper = 1) ## comparison


par(pin = c(5,3))
plot(x = 1:n, y = theta_est, type = "l", lty=2, ylim = c(0.7, 1.2), xlim = c(0, 10000))
polygon(x = c(1:n, rev(1:n)), y = c(parallelCI[,1], rev(parallelCI[,3])), lty = 4, border = "gray30", col = "gray90")
lines(x = 1:n, y = cumsum(h(z))/1:n, col = "red", lty = 3)

legend(x = 30000, y = 0.85, 
       legend = c(expression(hat(theta)), "CI", "New Experiment"), 
       border = "white", col = c("black", "blue", "red"), lty = c(2, 1, 3))




## An approximate but cheaper version of this basic Monte Carlo estimate of the variability is to bootstrap the originally obtained samples and from there estimate a 95% confidence band.
## bootstrap
M = 200L
ind = sample(x = 1:n, size = n*M, replace = TRUE) ## sample indices
boot_samples = matrix(data = h(x[ind]), nrow = n, ncol = M) ## matrix of samples
boot_est = apply(X = boot_samples, MARGIN = 2, FUN = cumsum) / 1:n ## matrix

bootCI = t(apply(X = boot_est, MARGIN = 1, FUN = quantile, c(0.025, 0.50, 0.975)))
summary(bootCI)




par(pin = c(5,3))
plot(x = 1:n, y = theta_est, type = "l", lty=2, ylim = c(0.7, 1.2), xlim = c(0, 10000))
polygon(x = c(1:n, rev(1:n)), y = c(bootCI[,1], rev(bootCI[,3])), lty = 4, border = "gray30", col = "gray80")
lines(x = 1:n, y = cumsum(h(z))/1:n, col = "red", lty = 3)

legend(x = 30000, y = 0.85, 
       legend = c(expression(hat(theta)), "CI", "New Experiment"), 
       border = "white", col = c("black", "blue", "red"), lty = c(2, 1, 3))



## Antithetic Variables  

n = 50000
u = runif(n)

h = function(x){exp(x)}

thetaHat = mean(h(u))
thetaHat
integrate(f = h, lower = 0, upper = 1)




## Antithetic Variables 


u1 = runif(n/2)
u2 = 1 - u1
h = function(x){exp(x)}

thetaHatAV = mean(c(h(u1), h(u2)))
thetaHatAV

cor(x = h(u1), y = h(u2)) ## correlation nearly -1











## Solution to Exercise 

n = 10^6
m = n/2

h <- function(x){ x/(2^x - 1) }

y = rnorm(n)
theta_MC = mean(h(y))

w = rnorm(m)
theta_AS = sum(h(w) + h(-1 * w)) / n

print(theta_MC)
print(theta_AS)




## Solution to Exercise 
## standard errors
rho = cor(h(w),h(-w))
se.a = (1+rho)*var(h(w))/n

print(rho)
print(se.a)



