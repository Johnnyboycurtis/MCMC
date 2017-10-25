

## The First Exercise  
set.seed(5678)
theta = 5 ## theta
rate = 1/theta ## R takes rate

t = 5 ## time cut off
N = 100 ## sample size of ex 1
M = 50 ## sample size of ex 2
y = rexp(n = N, rate = rate)
x = rexp(n = M, rate = rate)
x = sort(x)
E = as.integer(x > t)  ## 0 & 1

#N.ybar = sum(y)
ybar = mean(y)
Z = sum(E)
t = 5



## The First Exercise  

theta.j = 0.1
theta.jp1 = 0.5
for(i in 1:10){
  theta.j = theta.jp1
  p = (exp(-t/theta.j)/(1-exp(-t/theta.j)))
  theta.jp1 = (N*ybar + Z*( t + theta.j) + (M-Z)*(theta.j - t*p) ) / (N+M)
  print(theta.jp1)
}


## compare results
print(theta.jp1) ## EM theta estimate
mean(y) ## compare against MLE from observed data
mean(c(y, x)) ## compare against complete-data
## note, results will vary if you remove seed









## EM Normal Example  

set.seed(2345)
n = 100
mu = 4
sd = 1
x = rnorm(n, mu, sd) ## generate some data
c = 5 ## time cut off
w = x[x < c] ## obtain samples before time cut off
m = sum(x < c) ## number of observed samples
wbar = mean(w) ## observed mean
r = n - m ## difference in sample size


## EM Normal Example  
N = 200
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_old = mu_new
    mu_new = m*wbar/n + (r*mu_old/n) + 
      (r/n)*sd*(dnorm(c - mu_old))/(1 - pnorm(c - mu_old))  ## r/n instead of 1/n
    #print(mu_new)
}

print(tail(results))





## EM Normal Example  
plot(results, type = "l", main = "em estimates for mu", ylim = c(3.5, 4.5))
abline(h = mu, col = "red")
abline(h = wbar, col = "green", lty = 2)
abline(h = mean(x), col = "blue", lty = 3)



## Monte Carlo EM 

#A MC flavor of the EM algorithm

#1. Draw missing data sets $\mathbf{Z_1, Z_2, ..., Z_m} \sim f_{Z|X}(z | x, \theta_i)$ where each $\mathbf{Z_i}$ is a vector of all missing values needed to complete the observed data set $( \mathbf{X, Z} )$.

#2. Calculate $\bar{Q}(\theta | \theta_{i-1}, X, \mathbf{Z_1, ..., Z_m}) = \frac{1}{m} \sum_{i=1}^m Q(\theta | \theta_{i-1}, X, \mathbf{Z_i} )$

## EM Normal Example | Monte Carlo EM 

set.seed(2345)
n = 100
mu = 4
sd = 1
x = rnorm(n, mu, sd)
c = 5
w = x[x < c]
m = sum(x < c)
wbar = mean(w)
r = n - m


## EM Normal Example| Monte Carlo EM   

M = 10
N = 100
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_old = mu_new
    ## abs(N(0,1)) + mu_old + (c - mu_old) to *approximate*
    ## the truncated samples we need
    Z = matrix(data = (c - mu_old) + (mu_old +  abs(rnorm(n = r*M, mean = 0, sd = 1))), 
        nrow = r, ncol = M)
    mu_new = (m*wbar/n) + mean(colMeans(Z))*r/n
    M = M + 1
}


## EM Normal Example| Monte Carlo EM   


plot(results, type = "l", ylim = c(3.5, 4.5))
abline(h = mu, col = "red")
abline(h = wbar, col = "green", lty = 2)
abline(h = mean(x), col = "blue", lty = 3)






