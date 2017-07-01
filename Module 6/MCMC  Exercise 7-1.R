N = 10000

## For accept-reject, we need to find a value for M
## we can use `optimize` to find the maximum of our target density
maximum  = optimize(f = function(x){ dgamma(x = x, shape = 4.3, rate = 6.2)}, 
                    interval = c(0, 2), maximum = TRUE ) ## obtain maximum

M = maximum$objective
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), 
      from = 0, to = 2, add = TRUE, col = "blue",
      main = "Gamma(4.3, 6.2", xlab = "X", ylab = "Density")
abline(h = M, lty = 3, col = "red")

f = function(x){
  dgamma(x = x, shape = 4.3, rate = 6.2)
}

g = function(x){
  dgamma(x = x, shape = 4, rate = 7)
}


#N = 10
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

## see how samples from chain compare to Gamma density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), 
      from = 0, to = 2, col = "blue")






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



mean(X)

## see chain transitions
plot(X[500:1000], type = "o", main = "MCMC samples",
     xlab = "iterations", ylab = "X (samples obtained)")

## see how samples from chain compare to Gamma density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), 
      from = 0, to = 2, add = TRUE, col = "blue")



## Metropolis Hastings
## now compare results with Gamma(5,6)

N = 10000
X = numeric(N)
X[1] = rgamma(n = 1, shape = 4.3, rate = 6.2)
for(i in 1:N){
  Y = rgamma(n = 1, shape = 5, rate = 6)
  rho = (dgamma(x = Y, shape = 4.3, rate = 6.2) * dgamma(x = X[i], shape = 5, rate = 6)) / 
    (dgamma(x = X[i], shape = 4.3, rate = 6.2) * dgamma(x = Y, shape = 5, rate = 6)) 
  #X[i+1] = X[i] + (Y - X[i])*(runif(1) < rho) ## equivalent to if-else statement below
  if(runif(1) < rho){
    X[i+1] = Y
  } else{
    X[i+1] = X[i]
  }
  
}



mean(X)

## see chain transitions
plot(X[500:1000], type = "o", main = "MCMC samples",
     xlab = "iterations", ylab = "X (samples obtained)")

## see how samples from chain compare to Gamma density
hist(X, main = "Histogram of MCMC samples", prob = TRUE)
curve(expr = dgamma(x = x, shape = 4.3, rate = 6.2), 
      from = 0, to = 2, add = TRUE, col = "blue")

