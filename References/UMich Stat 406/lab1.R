x <- 3
print(x)
x <- 6 
y <- x
x1 <- 3; x2 <- 12; x3 <- 7; 
v <- c(x1, x2, x3) 
print(v) 
v <- c(3, 12, 7) 
print(v)
v <- seq(1,2,length=10)
print(v)
v <- seq(1,2,by=.1)
print(v)
v <- seq(0, 5, by=1) 
print(v)
w <- v^2 
print(w) 
w <- exp(v)
print(w)
print( round(w, 3) ) 
v <- c(4, 5, 10, 19, 23, 4, 16)  
print(v[6])
index <- c(1, 3, 4) 
print( v[index] ) 
index <- c(2,5,6) 
v[-index]
A <- matrix(0, nrow=3, ncol=3) 
print(A)
A[1,1] = 1 
A[1,2] = 4 
A[1,3] = 7
A[2,1] = 2 
A[2,2] = 5
A[2,3] = 8 
A[3,1] = 3 
A[3,2] = 6 
A[3,3] = 9
print(A) 
A[,2]
A[3,]
t(A)
v <- seq(0, 1, length=11) 
for(i in 1:11)
{
   print(v[i]) 
}
index <- c(1,2, 5, 8, 13, 22) 
for(i in index)
{
   print(i)
}
A = matrix(0, 5, 5) 
for(i in 1:5)
{ 
   for(j in 1:5)
   {
      A[i,j] = j^i 
   }
}
print(A)
if(x < 0) 
{
   print("x is negative") 
}  else
{
   print("x is positive") 
} 
if(x < 0) 
{
   print("x is negative") 
}  else if(x > 0)
{
   print("x is positive") 
}  else if(x == 0)
{
   print("x is 0")
}
k <- 0 

while(k <= 10) 
{

   print(k)
   k <- k + 1

}
m = 7
s <- seq(0, floor(m/2), by=1) 
k <- length(s)
bool <- FALSE
while( (bool == FALSE) & (k > 1) )
{

   bool <- ( (m/s[k]) == floor(m/s[k]) )
   k <- k - 1

}
sprintf("The largest factor of %i other than itself is %i", m, s[k+1])
n <- 10
a <- 1 
b <- 5 
u <- runif(n, min=a, max=b) 
n <- 10 
mu <- 0 
sigma <- 2
u <- rnorm(n, mean=mu, sd=sqrt(sigma))
n <- 10 
lambda <- 5
u <- rexp(n, rate=lambda) 
n <- 10 
k <- 20 
p <- .5 
u <- rbinom(n, size=k, prob=p) 
k <- 10
p <- .35
1-pbinom(3, size=k, prob=p)
prob <- 0 
for(x in 4:10) prob <- prob + choose(k, x)*(p^x)*((1-p)^(k-x))
nreps = 10000
X = rbinom(nreps, size=k,  prob=p) 
mean(X > 3)
lambda <- 5 
EX <- 1/lambda 
VarX <- 1/(lambda^2)
print(c(EX, VarX))
X <- rexp(10000, lambda) 
c(mean(X), var(X))  
mu <- 0 
sigma <- 1  
X <- rnorm(10000, mean=mu, sd=sqrt(sigma)) 
Y <- exp(X)/(1 + exp(X)) 
c( mean(Y), var(Y) ) 
theta <- seq(-5, 5, length=100) 
EY <- rep(0, 100)
VarY <- rep(0, 100)
for(i in 1:100)
{
 
      X <- rnorm(10000, mean=theta[i], sd=1) 
      
      Y <- exp(X)/(1 + exp(X)) 
      
      EY[i] <- mean(Y)
      VarY[i] <- var(Y)
      
}

plot( theta, EY,  xlab=expression(" " *theta), ylab="EY", main=expression("EY vs. " *theta), col=2 )
lines(theta,  EY, col=2)
 
plot( theta, VarY,  xlab=expression(" " *theta), ylab="VarY", main=expression("VarY vs. " *theta), col=3 )
lines(theta,  VarY, col=3)

