f <- function(x) .5-exp(-(x^2))
v <- seq(0,2,length=1000)
plt <- function(f,a0,b0)
{
   plot(v,f(v),type="l", ylim=c(-.5,.5), xlim=c(0,2),
   main=sprintf("Interval length = %f", (b0-a0)))
   abline(h=0)
   segments(a0,.1,a0,-.1,col=2,lty=1)
   points(a0,.1,col=2)
   points(a0,-.1,col=2)
   segments(b0,.1,b0,-.1,col=2,lty=1)
   points(b0,.1,col=2)
   points(b0,-.1,col=2)
   segments(mean(c(a0,b0)), .1, mean(c(a0,b0)), -.1, col=4, lty=1)
   points(mean(c(a0,b0)),.1,col=4)
   points(mean(c(a0,b0)),-.1,col=4)
}
iter1 <- function(f,I)
{
   m = mean(I)
   if( f(I[1])*f(m) < 0 ) return( c(I[1],m) ) else return( c(m, I[2]) ) 
}
I <- c(0,2)
f <- function(x) .5-exp(-(x^2))
a0 <- I[1]; b0 <- I[2]; 
plt(f, a0, b0) 
I = iter1(f,I)

bisect <- function(f, a, b, tol)
{
   I <- c(a,b) 
   L <- I[2]-I[1]
   while( L > tol ) 
   {
      m <- mean(I) 
      if( f(m)*f(I[1]) < 0 ) I = c(I[1],m) else I = c(m,I[2]) 
      L <- I[2]-I[1]
   }
   return(mean(I))
}
fn <- function(x) .5-exp(-(x^2))
bisect(fn,0,2,1e-6)
fn(.8325543)

X <- rnorm(100)  
fn <- function(t) sum(X-t) 
c( fn(-1), fn(1) ) 
bisect(fn, -1, 1, 1e-7)
mean(X)

f <- function(x) .5-exp(-(x^2)) 
df <- function(x) 2*x*exp(-(x^2))
v <- seq(0,2,length=1000)
plot(v,f(v),type="l",col=8)
abline(h=0)
x0 <- 1.5
segments(x0,0,x0,f(x0),col=2,lty=3)
slope <- df(x0) 
g <- function(x) f(x0) + slope*(x - x0)
v = seq(x0 - f(x0)/df(x0),x0,length=1000)
lines(v,g(v),lty=3,col=4)
x0 <- x0 - f(x0)/df(x0)

newton <- function(x0, f, df, d2f, tol=1e-4, pr=FALSE)  
{
   k <- 1
   fval <- f(x0) 
   grad <- df(x0) 
   hess <- d2f(x0) 
   xk_1 <- x0 
   
   cond1 <- sqrt( sum(grad^2) ) 
   cond2 <- Inf
   if( (cond1 < tol) ) return(x0) 
   
   while( (cond1 > tol) & (cond2 > tol) ) 
   {

      L <- 1
      bool <- TRUE 
      while(bool == TRUE) 
      {

         xk <- xk_1 - L * solve(hess) %*% grad
         if( f(xk) > fval ) 
         {
  
            bool = FALSE 
            grad <- df(xk) 
            fval <- f(xk)
            hess <- d2f(xk)  
         } else 
         {
           L = L/2
           if( abs(L) < 1e-20 ) return("Failed to find uphill step - try new start values")
         }

      }
      
      cond1 <- sqrt( sum(grad^2) ) 
      cond2 <- sqrt( sum( (xk-xk_1)^2 ))/(tol + sqrt(sum(xk^2)))
      k <- k + 1
      xk_1 <- xk 
     
   
   }

   if(pr == TRUE) print( sprintf("Took %i iterations", k) )
   return(xk)

}

## Ex 1 ## 
f <- function(x) exp(-(x^2)) 
df <- function(x) -2*x*f(x) 
d2f <- function(x) -2*f(x)-2*x*df(x) 
newton(2/3, f, df, d2f, 1e-7)

## Ex 2 ## 
X <- rnorm(100, mean=2, sd=sqrt(4))
n <- 100 
f <- function(t) 
{ 
    if( t[2] > 0) 
    { 
       return( sum( dnorm(X, mean=t[1], sd=sqrt(t[2]), log=TRUE) ) )
    } else 
    {
       return(-Inf) 
    }
}
df <- function(t) 
{ 
 mu <- t[1]; sig <- t[2]; 
 g <- rep(0,2) 
 g[1] <- (1/sig) * sum(X - mu) 
 g[2] <- (-n/(2*sig)) + sum( (X-mu)^2 )/(2*sig^2)
 return(g)
}
d2f <- function(t)
{
 mu <- t[1]; sig <- t[2]; 
 h <- matrix(0,2,2) 
 h[1,1] <- -n/sig
 h[2,2] <- (n/(2*sig^2)) - sum( (X-mu)^2 )/(sig^3)
 h[1,2] <- -sum( (X-mu) )/(sig^2)
 h[2,1] <- h[1,2]
 return(h)
}
newton( c(0,1), f, df, d2f)
c( mean(X), (n-1)*var(X)/n)

## Ex 3 ## 
n <- 100
X <- rgamma(n, shape=2, scale=3)
f <- function(t) 
{
    a <- t[1]; b <- t[2]; 
    if( (a > 0) & (b > 0) )
    {
       return( sum( dgamma(X, a, scale=b, log=TRUE) ) ) 
    } else
    {
       return(-Inf) 
    }
}
df <- function(t)
{
   g=c(0,0)
   a <- t[1]; b <- t[2]; 
   if( (a > 0) & (b > 0) )
   {
      g[1] <- -n*digamma(a) - n*log(b) + sum(log(X)) 
      g[2] <- -n*a/b + sum(X)/(b^2)
   } else
   {
     g <- c(-Inf, -Inf) 
   }
  return(g)
}
d2f <- function(t)
{
   h = matrix(0,2,2) 
   a <- t[1]; b <- t[2]; 
   if( (a>0) & (b>0) ) 
   {
      h[1,1] = -n*trigamma(a) 
      h[2,2] = n*a/(b^2) - 2*sum(X)/(b^3)
      h[1,2] = -n/b
      h[2,1] = -n/b
   } else
   {
     h <- matrix(-Inf, 2, 2)
   }
 return(h)
}
newton( c(1,5), f, df, d2f)
n <- 1e4
X = rgamma(n, shape=3, scale=7)
newton( c(1,5), f, df, d2f)

## Ex 4 ## 
n <- 100 
Beta <- c(.5, 1.2, -1.8, 3.6, 2.1) 
X <- cbind( rep(1,n), rnorm(n), rnorm(n), rnorm(n), rnorm(n) )
Y <- rep(0,n)
for(i in 1:n) Y[i] <- sum( Beta * X[i,] ) +  rnorm(1)

p <- ncol(X)  
L <- function(B) 
{
   like <- 0 
   for(i in 1:n)
   {
   
      mu_i <- sum( X[i,] * B ) 
      like <- like + dnorm(Y[i], mean=mu_i, sd=1, log=TRUE) 
      
   }
   
   return(like)
   
}

dL <- function(B) 
{
    grad <- rep(0, p) 
    for(i in 1:n)
    {
    
       mu_i <- sum( X[i,] * B )
       
       for(j in 1:p)
       {
       
          grad[j] <- grad[j] + X[i,j]*(Y[i] - mu_i)
          
       }
       
    }
    
    return(grad) 
    
}
d2L <- function(B) 
{

   H <- matrix(0, p, p) 
   for(i in 1:n)
   {
   
      for(j in 1:p)
      {
      
         for(k in 1:p)
         {
         
            H[j,k] <- H[j,k] - X[i,j]*X[i,k]
            
         }
         
      }
   
   }
   
   return(H)
   
}
newton(rep(0,p), L, dL, d2L)
lm(Y ~ X[,2] + X[,3] + X[,4] + X[,5])$coef
