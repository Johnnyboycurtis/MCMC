

u = runif(2000) ## generate a sample of unif(0,1) samples

hist(u, probability=TRUE, col = "gray", border = "white") ## plot histogram
## Q-Q plot for `runif` data against true theoretical distribution:
qqplot(x = qunif(ppoints(500)), y = u,  main = expression("Q-Q plot for Unif(0,1)"))
qqline(y = u, distribution = qunif, prob = c(0.1, 0.6), col = 2)

acf(x = u, main = "autocorrelation of samples")  ## autocorrelation function

## Generation of uniform samples 

require(graphics)
palette(value = rainbow(2))




qqplot(qunif(ppoints(500)),1-u, main = "U against 1 - U")
qqline(u, distribution = qunif, col = 2)

#It is important to generate samples that are not correlated with each other!



## Generation of uniform samples, pt. 2 
v = 10 * runif(1000, 0, 1) ## v ~ unif(0, 10)
hist(v, col = "gray", border = "white")
w = runif(1000, 0, 10) ## unif(0, 10)
hist(w, col = "gray", border = "white")






## Example: Binomial samples

#Another use of generating samples $U \sim U(0,1)$ is to generate Bernoulli. 
##Once we've generated $n$ Bernoulli samples, we can sum the Bernouilli samples to generate Binomial samples. Reference: [Notes](http://www.stat.ufl.edu/~abhisheksaha/sta4321/lect12.pdf)

# let's simulate 10,000 samples from bin(n=10, p = 0.3)
N = 10000 ## number of Unif and Binomial samples
n = 10 ## number of bernoulli trials for each bernoulli sample
p = 0.3 ## theoretical p
u =  runif(n = (n*N))
x = as.numeric(u <= p) ## convert bools to 1 and 0; Bernoulli samples
m = matrix(data = x , nrow = N, ncol = n)
Binom_samples = rowSums(m) ## binomial samples



## Example: Binomial samples, pt. 2

par(mfrow = c(1,2), mar = c(3, 4, 1, 2))
hist(Binom_samples, probability = TRUE, main = "Binom(10,0.3) from Unif(0,1)", 
     col = "blue", border = "white")
hist(rbinom(n = N, size = n, prob = p), probability = TRUE, main = "rbinom(10,0.3)",
     col = "gray", border = "white")
par(mfrow=c(1,1))






## Example: Exponential distribution 

N = 10^4
u = runif(N)

fInv = function(u){
  (-1/5) * log(u) ## or log(1-u)
}

outSamples = fInv(u)

## Example: Exponential distribution, pt. 2 
par(mfrow = c(1,2), mar = c(3, 4, 1, 2))
hist(outSamples, probability = TRUE, main = "Inv Trans", col = "gray", border = "white")
lines(x = ppoints(200), y = dexp(x = ppoints(200), rate = 5), 
      col = "blue")
hist(rexp(n = N, rate = 5), probability = TRUE, main = "rexp", col = "gray", border = "white", xlim = c(0, 1.5))
lines(x = ppoints(200), y = dexp(x = ppoints(200), rate = 5), 
      col = "blue")







## Example: Pareto Distribution, pt. 2 
set.seed(123)
n = 1000
U =runif(n)
a = 3
b = 2
X = b*(1-U)^(-1/a)
pareto = function(x){(a*(b^a)/x^(a+1))}

summary(X)


## Example: Pareto Distribution, pt. 3 
hist(X, probability = TRUE, breaks = 25, xlim =c(0, 20),
     col = "gray", border = "white",
     main = "Inverse Transform: Pareto(3,2)", xlab = "x")
curve(pareto(x), from = 0, to = 40, add = TRUE, col = "blue")






## Inverse Transform Discrete scenario, pt. 2 
n = 10000
u = runif(n)

fInv <- function(u){
  if(u <= 0.1) x <- 0
  if(u > 0.1 && u <= 0.3) x <- 1
  if(u > 0.3 && u <= 0.5) x <- 2
  if(u > 0.5 && u <= 0.7) x <- 3
  if(u > 0.7 && u <= 1) x <- 4
  return(x)
}

results = sapply(X = u, fInv)
table(results) / n



## Inverse Transform Discrete scenario, pt. 3 
z = table(results) / n
barplot(z, border = "white", col = "lightblue", xlab = "x", ylab = "probability", ylim = c(0, 0.4))








## Example: Beta(2,2) 

par(mfrow = c(1,2), mar = c(3, 4, 1, 2))
beta <- function(x) 6*x*(1-x)
M = 1.5
## plot 2
curve(expr = beta, from = 0, to = 1, 
      xlim = c(-0.5, 1.5), ylim = c(0,2),
      main = "Beta(2,2) Density with Unif(0,1)",
      xlab = "x", ylab = "density")

x = seq(from = -1, to = 2, by = 0.01)
Unif1 = function(x){ ifelse(x >= 0 & x <= 1, 1, 0) }
polygon(x, Unif1(x), lty = 9)

## plot 3
curve(expr = beta, from = 0, to = 1, 
      xlim = c(-0.2, 1.2), ylim = c(0,2),
      main = "Beta(2,2) with M*Unif(0,1)",
      xlab = "x", ylab = "density")

Unif2 = function(x){ ifelse(x >= 0 & x <= 1, 1*M, 0) }
polygon(x, Unif2(x), lty = 2)

#abline(h = 1.5, col = "red")
par(mfrow=c(1,1))







## Example: Beta(2,2), pt. 2
#Note that the target density $f$ has a maximum of 1.5, so we can set M = 1.5; 
##see: [Max of Beta(2,2)](http://www.wolframalpha.com/input/?i=max+6x(1+-+x))

## Accept-Reject
M = 1.5
X = rep(x = NA, 5) ## create a vector of length 5 of NAs
set.seed(123)
f <- function(x){ 6*x*(1 - x)} ## pdf of Beta(2,2)
g <- function(x){ 1 } ## pdf of Unif(0,1) is just 1

n = 10000

```


## Example: Beta(2,2), pt. 2

#First, we'll generate 5 samples to test the algorithm

for(i in 1:5){
  print(paste("Run: ", i))
  u = round(runif(1), 5)
  y = round(runif(1), 5)
  accept <- u <= f(y)/(M* g(y))
  print(paste("U: ", u, "and Y:", y, "and f/M*g: ",f(y)/(M* g(y))))
  print(paste("Accept? ", accept))
  if(accept){
    X[i] <- y
  }
}


## Example: Beta(2,2), pt. 3



X = rep(NA, n); M = 1.5
i = 0 ## index set to start at 0
while(sum(is.na(X))){
  U = runif(1); Y = runif(1)
  accept <- U <= f(Y)/(M*g(Y))
  if(accept){
    i = i+1 ## update the index
    X[i] <- Y
  }
}

round(summary(X), 4)
round(qbeta(p = c(0, 0.25, 0.5, 0.75, 1), 2, 2), 4)



## Example: Beta(2,2), pt. 4
par(mfrow = c(1,3), mar = c(3, 4, 1, 2))
beta <- function(x) 6*x*(1-x)

## plot 1
curve(expr = beta, from = 0, to = 1, 
      xlim = c(-0.5, 1.5), ylim = c(0,2),
      main = "Beta(2,2) Density with Unif(0,1)",
      xlab = "x", ylab = "density")

x = seq(from = -1, to = 2, by = 0.01)
Unif1 = function(x){ ifelse(x >= 0 & x <= 1, 1, 0) }
polygon(x, Unif1(x), lty = 9)

## plot 2
curve(expr = beta, from = 0, to = 1, 
      xlim = c(-0.2, 1.2), ylim = c(0,2),
      main = "Beta(2,2) with M*Unif(0,1)",
      xlab = "x", ylab = "density")

Unif2 = function(x){ ifelse(x >= 0 & x <= 1, 1*M, 0) }
polygon(x, Unif2(x), lty = 2)

## plot 3
hist(X, xlab = "X", main = "Beta(2,2) from Accept-Reject algorithm", 
     probability = TRUE, col = "gray", border = "white")
curve(expr = beta, from = 0, to = 1, add = TRUE)







## Example: Beta(2,2)

N = 1000; U = runif(N); Y = runif(N); M = 1.5

f <- function(x){ 6*x*(1 - x)} ## pdf of Beta(2,2)
g <- function(x){ 1 } ## pdf of Unif(0,1) is just 1

accept <- U*M < f(Y)/(g(Y))

mean(accept) ## acceptance rate
print(1/M) ## probability of acceptance



##Exercise: Change the value of M to see how the performance of the algorithm changes

## Example: Beta(2,2)

par(mar = c(3, 4, 1, 2))
palette(value = rainbow(6))
plot(Y, U*M, col = as.numeric(accept)+3, 
     xlim = c(-0.2, 1.2), ylim = c(0,2),
     main = "Accept-Reject with M = 1.5")

curve(expr = beta, from = 0, to = 1, 
      xlim = c(-0.5, 1.5), ylim = c(0,2),
      main = "Beta(2,2) with M*Unif(0,1)",
      xlab = "x", ylab = "density", add = TRUE,
      lwd = 4)

Unif2 = function(x){ ifelse(x >= 0 & x <= 1, 1*M, 0) }
polygon(x, Unif2(x), lty = 2, lwd = 2)




#It should be noted that the probability of acceptance is given by $\frac{1}{M}$. So, in order to make an efficient accept-reject algorithm, we should set $M$ to be as high as needed, but no larger! As $M$ increases, the probability of acceptance decreases, and this results in an increase in draws where we do not obtain samples from our target distribution $f$. This increases the computational cost. 








N = 10000
Z = rnorm(n = 2*N, mean = 0, sd = 1)
Z = matrix(data = Z, nrow = N, ncol = 2)

transfromation <- function(vec){
  R = sqrt(sum(vec^2))
  #R = sqrt(vec[1]^2 + vec[2]^2)
  return(R)
}

R_Out = apply(X = Z, MARGIN = 1, FUN = transfromation)

summary(R_Out)
sqrt(pi/2) ## theoretical mean





## Example: Generate Rayleigh samples
hist(R_Out, col = "gray", border = "white")
## compare with Rayleigh {VGAM}	








# Appendix

## A Solution to Beta transformation


par(mar = c(3, 4, 1, 2))
a = 3                                                    
b = 10

A = matrix(rexp(10000*a, 1), ncol = a)
B = matrix(rexp(10000*b, 1), ncol = b)
numerator = rowSums(A)
denom = rowSums(A) + rowSums(B)
betasamples = numerator/denom

hist(betasamples, probability=TRUE, col = "gray", border = "white")
curve(dbeta(x, a, b), from=0, to=1, add=TRUE, lty=3, lwd = 2)









n = 10000
x = rexp(n=n, rate = 1)
mat = matrix(data = x, ncol = 4)
chisqResults = rowSums(mat)
par(mfrow = c(1, 3))
hist(chisqResults)
hist(rchisq(n, df = 2*4))
qqplot(x = chisqResults, y = rchisq(n, 8))



b = 4
a = 2
n = 10000
x = rexp(n, 1)
mat = matrix(x, ncol = a)
gammaResults = b*rowSums(mat)
qqplot(x = gammaResults, y = rgamma(n,a,b))

