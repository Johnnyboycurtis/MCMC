## standard normal mc em example 6.2

set.seed(345)

N = 100
C = 5
mu = 4
sd = 1

samples = sort(rnorm(N, mu, sd))

M = sum(samples > C)
wbar = mean(samples[samples <= C])


mu_kplus1 =  function(n, m, wbar, mu_k, c){
    #out = (m*wbar/n) + (n-m)*mu_k/n + (1/n)*(dnorm(c-mu_k) / (1 - pnorm(c - mu_k)))
    out = (m*wbar/n) + (n-m)*mu_k/n + (1/n)*(dnorm(c-mu_k) / (1 - pnorm(c - mu_k)))
    return(out)

}


mu_kp1 = 1

for(i in 1:12){
    mu_k = mu_kp1
    mu_kp1 <- mu_kplus1(N, M, wbar, mu_k, C)
    print(mu_kp1)
}

mu_k = 0
mu_kp1 = 1
results = numeric(1000)
for(i in 1:1000){
    mu_k = mu_kp1
    mu_kp1 <- mu_kplus1(N, M, wbar, mu_k, C)
    results[i] <- mu_kp1
}


plot(results, type = "l")













set.seed(123)
n = 100
mu = 4
sd = 1
x = rnorm(n, mu, sd)
c = 5

w = x[x < c]
r = sum(x < c) ## number of observed
wbar = mean(w) ## mean of observed
m = sum(x>c) ## number of missing

M = 5
N = 10
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_old = mu_new
    Y = matrix(data = rep(y, M), nrow = r, ncol = M)
    Z = c + matrix(data = rnorm(n = m*M, mean = mu_old, sd = 1), 
        nrow = m, ncol = M)
    X = rbind(Y,Z)
    mu_new = r*wbar/n + ( m*mean(colMeans(Z))/n ) 
    #M = M + 1
}


print(tail(results))
plot(results, type = "l", main = "em estimates for mu")






##############################################################################



n = 100
mu = 4
sd = 1
x = rnorm(n, mu, sd)
c = 5
w = x[x < c]
m = sum(x < c)
wbar = mean(w)
r = n - m

N = 10
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_k = mu_new
    mu_new = (m*wbar/n) + (n-m)*mu_k/n + (1/n)*(dnorm(c-mu_k) / (1 - pnorm(c - mu_k)))
    print(mu_new)
}

plot(results, type = "l", ylim = c(3.5, 4.5))
abline(h = mu, col = "red")
abline(h = wbar, col = "green", lty = 2)
abline(h = mean(x), col = "blue", lty = 3)






## working example as well
r = n - m
M = 10
N = 1000
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_old = mu_new
    Y = matrix(data = rep(w, M), nrow = m, ncol = M)
    Z = matrix(data = (c - mu_old) + rnorm(n = r*M, mean = mu_old, sd = 1), 
        nrow = r, ncol = M)
    X = rbind(Y,Z)
    mu_new = mean(colSums(X))/n
    #print(colMeans(X))
    #print(mu_new)
    #mu_new = r*wbar/n + ( m*mean(colMeans(Z))/n ) 
    M = M + 1
}

plot(results, type = "l")
abline(h = mu, col = "red")




## ideal example
n = 100
mu = 4
sd = 1
x = rnorm(n, mu, sd)
c = 5
w = x[x < c]
m = sum(x < c)
wbar = mean(w)
r = n - m

M = 10
N = 2000
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_old = mu_new
    Z = matrix(data = (c - mu_old) + rnorm(n = r*M, mean = mu_old, sd = 1), 
        nrow = r, ncol = M)
    mu_new = (m*wbar/n) + mean(colMeans(Z))*r/n
    M = M + 1
}

plot(results, type = "l", ylim = c(3.5, 4.5))
abline(h = mu, col = "red")
abline(h = wbar, col = "green", lty = 2)
abline(h = mean(x), col = "blue", lty = 3)






## alternatively, abs(rnorm(n, 0, 1)) simulations to generate right truncated sims

M = 10
N = 200
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_old = mu_new
    Z = matrix(data = (c - mu_old) + (mu_old +  abs(rnorm(n = r*M, mean = 0, sd = 1))), 
        nrow = r, ncol = M)
    mu_new = (m*wbar/n) + mean(colMeans(Z))*r/n
    M = M + 1
}

plot(results, type = "l", ylim = c(3.5, 4.5))
abline(h = mu, col = "red")
abline(h = wbar, col = "green", lty = 2)
abline(h = mean(x), col = "blue", lty = 3)


