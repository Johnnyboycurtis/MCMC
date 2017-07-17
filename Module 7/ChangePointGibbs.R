
data_array = c(4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6,
                    3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
                    2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0,
                    1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
                    0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2,
                    3, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1)

disasters = ts(data_array, freq=1, start = 1851)

plot(disasters, type = "o", lty = 2)

plot(disasters, type = "h", lty = 6, lwd = 2)


y = data_array



# Gibbs sampler for the coal mining change point
# initialization
n <- length(y)
 #length of the data
m <- 10^4
 #length of the chain
mu <- numeric(m)
lambda <- numeric(m)
k <- numeric(m)
L <- numeric(n)
k[1] <- sample(1:n, 1)
mu[1] <- 1
lambda[1] <- 1
b1 <- 1
b2 <- 1




# run the Gibbs sampler
for (i in 2:m){
    kt <- k[i-1]
    #generate mu
    r <- .5 + sum(y[1:kt])
    mu[i] <- rgamma(1, shape = r, rate = kt + b1)
    #generate lambda
    if (kt + 1 > n) r <- .5 + sum(y) else
    r <- .5 + sum(y[(kt+1):n])
    lambda[i] <- rgamma(1, shape = r, rate = n - kt + b2)
    #generate b1 and b2
    b1 <- rgamma(1, shape = .5, rate = mu[i]+1)
    b2 <- rgamma(1, shape = .5, rate = lambda[i]+1)


    for (j in 1:n) {
        L[j] <- exp((lambda[i] - mu[i]) * j) * 
                (mu[i] / lambda[i])^sum(y[1:j])
    }

    L <- L / sum(L)
    #generate k from discrete distribution L on 1:n
    k[i] <- sample(1:n, prob=L, size=1)
}







b <- 201
j <- k[b:m]

print(mean(k[b:m]))
#[1] 39.935

print(mean(lambda[b:m]))
#[1] 0.9341033

print(mean(mu[b:m]))
#[1] 3.108575



## Code for Figure 9.12 on page 276
# histograms from the Gibbs sampler output
par(mfrow=c(2,3))
labelk <- "changepoint"
label1 <- paste("mu", round(mean(mu[b:m]), 1))
label2 <- paste("lambda", round(mean(lambda[b:m]), 1))

hist(mu[b:m], main="", xlab=label1,
breaks = "scott", prob=TRUE) #mu posterior

hist(lambda[b:m], main="", xlab=label2,
breaks = "scott", prob=TRUE) #lambda posterior

hist(j, breaks=min(j):max(j), prob=TRUE, main="",
xlab = labelk)

par(mfcol=c(1,1), ask=FALSE) #restore display