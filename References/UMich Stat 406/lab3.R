###### Section 1 #######

x <- seq(-3, 3, length=1000)
phi <- rep(1/sqrt(2*pi), 1000)
for(i in 1:1000)
{

  phi[i] <- phi[i]*exp(-(x[i]^2)/2 ) 

}
f <- function(x) 
{
  return( exp(-(x^2)/2) ) 
}

x <- seq(-3, 3, length=1000) 
phi <- rep(1/sqrt(2*pi), 1000)

for(i in 1:1000)
{
 
  phi[i] <- phi[i]*f(x[i])

}
# global variable 
x <- 5

# x is also the input name
f <- function(x) 
{
  print(x) 
}
f(2)
f <- function(y)
{

  print(x)
  x <- y  
  print(x)

}
f(2)
f <- function(z,y) 
{
  return( x + y ) 
}
f(2,3) 
x <- 10 
f <- function(y) 
{
 print(x) 
 x = x + y 
 print(x) 
}
f(5) 
print(x) 
get.root <- function(a, b, c)
{

  d <- b^2 - 4*a*c 

  if( d < 0 ) 
  {

    roots <- c()

  } else
  {

    roots <- c( (-b - sqrt(d))/(2*a), (-b + sqrt(d))/(2*a) ) 
  
  }

  return(roots)

}
get.root(1, 1, -2) 
get.root(1, 0, 1)

####### Section 2 ########

A = matrix( rbinom(2000, 50, .6), 200, 10) 
f <- function(x) sort(x)[3] 
OUT <- rep(0, 200)  
for(i in 1:200)
{

   OUT[i] <- f( A[i,] ) 

}
 
OUT <- apply(A, 1, f) 
f <- function(x) x[1]^2 + x[2]^2 
x1 <- seq(0,1,by=.01)
x2 <- seq(1,2,by=.01) 
q <- NULL
k <- 0 
for(i in 1:101)
{
 for(j in 1:101)
 {
  k=k+1
  q[k] <- f( c(x1[i],x2[j]) ) 
 }
}
Q <- t( matrix(q, 101, 101) )
n <- NULL
for(i in 1:101) n <- c(n, rep((i-1)/100, 101)) 
v <- cbind(n, seq(1,2,length=101))
z <- apply(v,1,f)
Z <- t( matrix(z, 101, 101) )

####### Section 3 ########

rproc <- function(n, p.i, p.d) 
{

   P <- matrix(0, 3, 3) 
   P[1,] <- c(1-p.i, p.i, 0) 
   P[2,] <- c(.5*(1-p.d), .5*(1-p.d), p.d)
   P[3,] <- c(0, 0, 1) 

   X <- list() 

   for(j in 1:n)
   {

      x <- c(1)
      k <- 1 
      
      while( x[k] != 3 ) 
      {

          u <- runif(1) 
          p <- P[x[k],]
          if( u < p[1] )
          {

             x <- c(x,1) 

          } else if( u < sum(p[1:2]) ) 
          {
          
             x <- c(x,2) 

          } else 
          {
           
             x <- c(x,3)
             
          }
          
          k <- k + 1 
          
       }
       
       X[[j]] <- x	     
    }
   return(X) 
}
char.proc <- function(Q)
{ 
   n <- Q[1]; p.i <- Q[2]; p.d <- Q[3]; T <- Q[4]; 
   Y <- rproc(n, p.i, p.d) 
   L <- rep(0, n)
   for(i in 1:n)
   {   
      L[i] <- length( Y[[i]] )  
   }
   return( c( mean(L <= T), mean(L) ) ) 
}
pi <- seq(.1, .9, by=.2)
pd <- seq(.1, .9, by=.2)  
G <- NULL
for(i in 1:5) G <- c(G, rep(pi[i], 5)) 
G <- cbind(G, pd)
G <- cbind(rep(200,25), G, rep(5, 25) ) 
OUT <- apply(G, 1, char.proc)
ii1 <- seq(1, 25, by=5) 
ii2 <- seq(2, 25, by=5)
ii3 <- seq(3, 25, by=5)
ii4 <- seq(4, 25, by=5)
ii5 <- seq(5, 25, by=5)
plot( G[ii1,2], OUT[1,ii1], col=2, ylim=c(0,1), xlab="p.i", 
   ylab="Probability of death by time 5", axes=F, 
   main="Prob. of death by time 5 vs. p.i for each p.d")
axis(1, seq(.1, .9, length=5), seq(.1, .9, length=5))
axis(2, seq(0, 1, length=6), seq(0, 1, length=6))
box()
lines( G[ii1,2], OUT[1,ii1], col=2) 
lines( G[ii2,2], OUT[1,ii2], col=3)
points( G[ii2,2], OUT[1,ii2], col=3)
points( G[ii3,2], OUT[1,ii3], col=4)
lines( G[ii3,2], OUT[1,ii3], col=4)
points( G[ii4,2], OUT[1,ii4], col=5)
lines( G[ii4,2], OUT[1,ii4], col=5)
points( G[ii5,2], OUT[1,ii5], col=6) 
lines( G[ii5,2], OUT[1,ii5], col=6) 
plot( G[ii1,2], OUT[2,ii1], col=2, ylim=c(0,80), xlab="p.i", 
   ylab="Expected time until death", axes=F, 
   main="Expected time until death vs. p.i for each p.d")
axis(1, seq(.1, .9, length=5), seq(.1, .9, length=5))
axis(2, seq(0, 80, length=5), seq(0, 80, length=5))
box()
lines( G[ii1,2], OUT[2,ii1], col=2) 
lines( G[ii2,2], OUT[2,ii2], col=3)
points( G[ii2,2], OUT[2,ii2], col=3)
points( G[ii3,2], OUT[2,ii3], col=4)
lines( G[ii3,2], OUT[2,ii3], col=4)
points( G[ii4,2], OUT[2,ii4], col=5)
lines( G[ii4,2], OUT[2,ii4], col=5)
points( G[ii5,2], OUT[2,ii5], col=6) 
lines( G[ii5,2], OUT[2,ii5], col=6) 

######## Section 4 #########

bet.round <- function(I) 
{

   p <- I[1]; D <- I[2]; 
   startval <- D
   bank <- NULL
   k <- 0 
   while( (D > 0) & (D < (2*startval)) ) 
   {
      ind <- rbinom(1, 1, p) 
      if( ind == 1) 
      {
      
         D = D + 1 
        
      } else 
      {
      
         D = D - 1
         
      }
      
      k <- k + 1
      bank = c(bank, D) 
      
   }
   
   # return(c(k, (D==0)))
   
   return(bank) 
   
}
plot( bet.round( c(.4, 50) ), xlab="Time", ylab="Bankroll", 
   main="Random bettor's trajectory", type="l")
bet.round <- function(I) 
{
   p <- I[1]; D <- I[2]; 
   startval <- D
   bank <- NULL
   k <- 0 
   while( (D > 0) & (D < (2*startval)) ) 
   {
      ind <- rbinom(1, 1, p) 
      if( ind == 1) 
      {     
         D = D + 1       
      } else 
      {     
         D = D - 1        
      }     
      k <- k + 1
      bank = c(bank, D)     
   }  
   return(c(k, (D==0)))
   #return(bank) 
}
meantime <- function(J)
{
 
   p <- J[1]; D <- J[2]; nreps <- J[3]; 
   m <- matrix(0, nreps, 2)
   for(j in 1:nreps) m[j,] <- bet.round( c(p, D) ) 
   return( c(mean(m[,1]), mean(m[,2]) ) )
   
}
probs <- as.matrix( seq(.4, .6, length=10), 1, 10) 
inputs <- cbind( probs, rep(10, 10), rep(1000, 10) ) 
times <- apply(inputs, 1, meantime)
plot(probs, times[2,], xlab="Probability of winning a single round", 
   ylab="Probability of ruin", 
   main="Probability of ruin vs. Probability of winning a single round", col=3)
lines(probs, times[2,], col=3)
plot(probs, times[1,], xlab="Probability of winning a single round", 
   ylab="Expected time until ending", 
   main="Probability of ruin vs. Expected time until you quit playing", col=4)
lines(probs, times[1,], col=4)

######## Section 5 ##########

single.tau <- function(k) 
{
   X <- c(0) 
   i <- 0 
   while( sum(abs(X)) < k ) 
   {
      i <- i + 1 
      X<- c(X, rnorm(1)) 
   } 
   return(i)
}
rtau <- function(nreps, k)
{
   B <- as.matrix( rep(k, nreps), 1, nreps ) 
   V <- apply(B, 1, single.tau)
   return(V)
}
K <- seq(.5, 5, by=.5) 
Q <- matrix(0, 10, 2)
for(j in 1:10)
{

   TAU <- rtau(1000, K[j])
   Q[j,] <- c( mean(TAU), sd(TAU) ) 
   
}

par(mfrow=c(2,1))
plot(K, Q[,1], xlab="The cut off point, k", 
   ylab="Mean of Tau", main="Mean of tau vs. k", col=2)
lines(K, Q[,1], col=2)

plot(K, Q[,2], xlab="The cut off point, k", 
   ylab="S.D. of Tau", main="S.D. of tau vs. k", col=3)
lines(K, Q[,2], col=3)
