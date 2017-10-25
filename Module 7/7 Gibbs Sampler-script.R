## Bivariate Normal Example 
library(MASS, quietly = TRUE) ## for kernel density estimation
library(RColorBrewer, quietly = TRUE) ## some pretty colors
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

N = 100000
rho = 0.9
mu_x = 1
mu_y = 2
sd_x = 1.2
sd_y = 0.75

s1 = sqrt(1 - rho^2) * sd_x
s2 = sqrt(1 - rho^2) * sd_y

MVN = matrix(data = NA, nrow = N, ncol = 2, 
           dimnames = list(NULL, c("X", "Y")))




## Bivariate Normal Example 
MVN[1, ] = c(mu_x, mu_y)
Y = MVN[1, 2] ## get Y vals
for(i in 1:(N-1)){
  mx = mu_x + rho * (Y - mu_y) * sd_x/sd_y
  X = rnorm(n = 1, mx, s1)
  MVN[i+1, 1] = X
  my = mu_y + rho * (X - mu_x) * sd_y/sd_x
  Y = rnorm(n = 1, mean = my, sd = s2)
  MVN[i+1, 2] = Y
}


## means
colMeans(MVN)

## correlation
cor(MVN)





## Bivariate Normal Example 
z = kde2d(x = MVN[,1], y = MVN[,2], n = 50) ## kernel density estimation
plot(MVN, xlab="X", ylab="Y", pch = 18, col = 8, cex=.4, main = "MVN Samples")
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE) ## add contours to plot






## Beta-Binomial revisited 
## Beta-Binomial revisited 
n = 16
a = 2
b = 4
X = numeric(N)
thetas = numeric(N)

X[1] = runif(1) ## initial values
thetas[1] = runif(1) ## theta values

for(i in 1:N){
  X[i+1] = rbinom(n = 1, size = n, prob = thetas[i])
  thetas[i+1] = rbeta(1, a + X[i+1], n - X[i+1]+b)
}

quantile(X)
quantile(thetas)




## Beta-Binomial revisited 



par(mfrow = c(1,2))
hist(X, main = "Marginal Dist: X", probability = TRUE, border = "white", col = "gray")
hist(thetas, main = paste("Marginal Dist:", expression(theta)), probability = TRUE, border = "white", col = "gray")

## Bayesian Change-Point Analysis 

data_vector = c(4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 
                3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
                2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0, 
                1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
                0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 
                3, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1)
n = length(data_vector)
disasters = ts(data_vector, freq=1, start = 1851)
plot(disasters, type = "h", lty = 6, lwd = 2, main = "Coal Mining Disasters", xlab = "Year", ylab = "Disasters")





## Bayesian Change-Point Analysis 
## Implementing Gibbs sampling  
set.seed(123)

y = data_vector
# Gibbs sampler for the coal mining change point
# initialization
n <- length(y) #length of the data
m <- 10^4 #length of the chain

## vectors to hold data
theta <- numeric(m) 
lambda <- numeric(m) 
k <- numeric(m)
L <- numeric(n)

## initial values
k[1] <- sample(1:n, 1) ## change-points
theta[1] <- 1
lambda[1] <- 1
a = 0.5
b1 <- 1
b2 <- 1


# run the Gibbs sampler
for (t in 2:m){
    kt <- k[t-1]
    #generate theta
    r <- a + sum(y[1:kt])
    theta[t] <- rgamma(1, shape = r, rate = kt + b1)
    #generate lambda
    if (kt + 1 > n){ 
      r <- a + sum(y) 
      }else{
        r <- a + sum(y[(kt+1):n])
      }
    lambda[t] <- rgamma(1, shape = r, rate = n - kt + b2)
    #generate b1 and b2
    b1 <- rgamma(1, shape = a, rate = theta[t]+1)
    b2 <- rgamma(1, shape = a, rate = lambda[t]+1)

    for (j in 1:n) {
        L[j] <- exp((lambda[t] - theta[t]) * j) * 
                (theta[t] / lambda[t])^sum(y[1:j])
    }
    L <- L / sum(L)
    #generate k from discrete distribution L on 1:n
    k[t] <- sample(1:n, prob=L, size=1)
}


#Set a burn-in of 1000 samples. We will use burn in to toss out "poor" samples from out Markov chain. Arguments for and against burn-in vary. Statisticians, Andrew Gelman ([Burn-in Man](http://andrewgelman.com/2008/09/06/burnin_man/)) and Charlie Geyer ([Burn-In](http://users.stat.umn.edu/~geyer/mcmc/burn.html)) provide some commentary on burn-in.
## set burn in
burn_in <- 1000
## will toss out first 1000 samples

K <- k[burn_in:m]

## mean
print(mean(K))
#[1] 39.935

## mode
print(which.max(table(K)))



## Implementing Gibbs sampling  


plot(K, type="l", col="gray", main = "Change-Point",
     xlab="Iteration", ylab = "K")
lines(1:length(K), cumsum(K) / (1:length(K)))




## Implementing Gibbs sampling  


print(mean(lambda[burn_in:m])) #[1] 0.9341033
print(mean(theta[burn_in:m])) #[1] 3.108575




## Code for Figure 9.12 on page 276
# histograms from the Gibbs sampler output
par(mfrow=c(1,3))

hist(theta[burn_in:m], main="", xlab = expression(theta), 
     col="gray", border="white",
     breaks = "scott", prob=TRUE) #theta posterior

hist(lambda[burn_in:m], main="", xlab = expression(lambda),
     col="gray", border="white",
     breaks = "scott", prob=TRUE) #lambda posterior

hist(K, breaks=min(K):max(K), prob=TRUE, main="",
     col="gray", border="white",
     xlab = "k")


