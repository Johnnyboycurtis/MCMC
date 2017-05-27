rX <- function(n, p) 
{
   U <- runif(n) 
   X <- rep(0,n) 
   w1 <- which(U <= p[1])
   X[w1] <- 1 
   w2 <- which( (U > p[1]) & (U < sum(p[1:2])) )
   X[w2] <- 2 
   w3 <- which( U > sum(p[1:2]) )
   X[w3] <- 3 
   return(X)
}
X = rX(10000, c(.4, .25, .35))
mean(X == 1) 
mean(X == 2) 
mean(X == 3) 
pois.cdf <- function(x, L) 
{
 
  # grid of values at which to calculate the mass function 
  v <- seq(0, x, by=1)
  
  # return CDF value
  return( sum( exp(-L)*(L^v)/factorial(v) ) ) 

}
r.pois <- function(n, L) 
{

   U <- runif(n) 
   X <- rep(0,n) 
   for(i in 1:n) 
   {
     if(U[i] < pois.cdf(0,L)) 
     {

        X[i] <- 0

     } else
     { 
        B = FALSE 
        I = 0 
        while(B == FALSE) 
        {
           int <- c( pois.cdf(I, L), pois.cdf(I+1,L) ) 
           if( (U[i] > int[1]) & (U[i] < int[2]) ) 
           {
              X[i] <- I+1
              B = TRUE

           } else 
           { 
             I=I+1
           }               


       }

    }


  }
  
  return(X)

} 
V = r.pois(1000, 4) 
Mf <- c(0:15)
for(i in 0:15) Mf[i+1] <- mean(V==i) 
b <- c(0:15)
plot(b, Mf, xlab="x", ylab="p(x)", main="Empirical mass function", col=2)
lines(b, Mf, col=2)
points(b, dpois(b,4), col=4)
lines(b, dpois(b,4), col=4)
r.harmonic <- function(n)
{
  
   # generate uniforms
   U <- runif(n)
   
   # return F^-1(U)
   return( U/(1-U) )

}
X <- r.harmonic(1000)
v <- seq(0, 20, length=1000) 
emp.cdf <- rep(0, 1000) 
for(i in 1:1000) emp.cdf[i] <- mean(X <= v[i])
true.cdf <- v/(1+v) 
plot(v, emp.cdf, xlab="X", ylab="F(X)", main="Empirical vs True CDF", col=2,type="l")
lines(v,true.cdf,col=4) 
r.skewlog <- function(n, Theta)
{

    U <- runif(n) 
    return( -log( (1/U)^(1/Theta) - 1) )

}

X <- r.skewlog(1000, 4) 
v <- seq(-10, 10, length=1000) 
emp.cdf <- rep(0, 1000) 
for(i in 1:1000) emp.cdf[i] <- mean(X <= v[i]) 
true.cdf <- (1 + exp(-v))^(-4) 
plot(v, emp.cdf, xlab="X", ylab="F(X)", main="Empirical vs True CDF", col=2,type="l")
lines(v,true.cdf,col=4) 
X <- r.skewlog(1000, 3) 
Y <- log(1 + exp(-X)) 
Lambda <- seq(.1, 5, length=1000) 
M <- rep(0, 1000) 
for(j in 1:1000)
{
   G <- function(x) Lambda[j]*( (1+x)^2 )*exp(-Lambda[j]*x)
   if( Lambda[j] < 2 )
   {
      M[j] <- G(-1 + 2/Lambda[j])
   } else
   {
      M[j] <- G(0)
   }
}
plot(Lambda, M, type="l", xlab="Lambda", ylab="M", main="M as a function of Lambda") 
Lambda = 1.5
N = 1000
p1 <- function(x) Lambda*exp(-x*Lambda) 
p2 <- function(x) (1 + x)^(-2)
r.trial <- function(n) 
{
   U <- runif(n) 
   return( U/(1-U) ) 
}
if(Lambda < 2)
{ 
   M = p1(-1 + 2/Lambda)/p2(-1 + 2/Lambda) 
} else
{
   M = p1(0)/p2(0)
}
#M = 100
X <- rep(0, N) 
niter = 0 

for(j in 1:N)
{
   B <- FALSE 

   while(B == FALSE)
   {  
      niter = niter + 1
      V <- r.trial(1) 
      R <- p1(V)/(M*p2(V)) 
      U <- runif(1) 
      if(U < R) 
      {

         X[j] <- V   
         B <- TRUE 

      }

   }

}

M
mean(X)
var(X)
niter
M
mean(X)
var(X)
niter
