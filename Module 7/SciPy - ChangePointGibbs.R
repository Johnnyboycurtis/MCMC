
data_array = c(4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6,
                    3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
                    2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0,
                    1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
                    0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2,
                    3, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1)

disasters = ts(data_array, freq=1, start = 1851)

#plot(disasters, type = "o", lty = 2)

#plot(disasters, type = "h", lty = 6, lwd = 2)


y = data_array




alpha = 0.6
beta = 10

# Specify number of iterations
N = 6000
B = 1000

# Initialize trace of samples
mu = numeric(N+1)
lambda2 = numeric(N+1)
tau = numeric(N+1)


## initialize values
mu[1] = 6
lambda2[1] = 2
tau[1] = 50


n = length(disasters)

DGamma = function(lambda, a, b){
  lambda**(a-1) * exp(-b*lambda)
}





y = disasters
# Sample from conditionals
for(i in 1:N){
    kt = tau[i]
  # Sample early mean
  if (kt + 1 > n){
      r <- alpha + sum(y)
   }else{
    r = sum(y[1:kt]) + alpha ## disasters_array[tau[i]:].sum()
   }
  mu[i+1] = rgamma(n = 1, shape = r, rate = kt + 1)

  # Sample late mean
  lambda2[i+1] = rgamma(n = 1, shape = r, rate = (n - kt + 1))

  # Sample changepoint: first calculate probabilities (conditional)
  p = numeric(n)

#   for(j in 1:n){
#       G1 = DGamma(mu[i+1], sum(y[1:j]) + alpha, j + beta)
#       G2 = DGamma(lambda2[i+1], sum(y[(j):111]) + alpha, n - j + beta)
      
#       p[j] <- exp(G1) + exp(G2)
#   }

#   if(i < 4){ 
#       #print(G1)
#       #print(G2)
#       print(p) 
#       print(length(p))
#       }
#   # # ... then draw sample
#   tau[i+1] = sample(x = 1:111, size = 1, prob = p/sum(p))
    for (j in 1:n) {
            p[j] <- exp((lambda2[i+1] - mu[i]) * j) * 
                    (mu[i] / lambda2[i+1])^(sum(y[1:j]))
        print(p[j])
        }
    if(sum(is.na(p))) break
    p <- p / sum(p)
    #generate k from discrete distribution L on 1:n
    tau[i+1] <- sample(1:n, prob=p, size=1)
    #break
}




tauPost = tau[-(1:B)]
lambdaPost1 = mu[-(1:B)]
lambdaPost2 = lambda2[-(1:B)]

#plot(density(tauPost))

#plot(density(lambdaPost1))

#plot(density(lambdaPost2))

