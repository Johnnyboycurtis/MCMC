
set.seed(123)
n = 100
mu = 4
sd = 3
x = rnorm(n, mu, sd)
c = 7

w = x[x<c]
wbar = mean(w)
m = sum(x>c)


N = 200
mu_new = wbar
results = numeric(N)
for(i in 1:N){
    results[i] = mu_new
    mu_old = mu_new
    mu_new = m*wbar/n + (n - m)*mu_old/n + (1/n)*sd*(dnorm(c - mu_old))/(1 - pnorm(c - mu_old)) 
    #print(mu_new)
}


print(tail(results))
plot(results, type = "l", main = "em estimates for mu", ylim = c(3,5))





