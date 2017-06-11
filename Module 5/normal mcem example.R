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
