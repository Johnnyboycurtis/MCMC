k = 1000 
n = 50 
THETA <- seq(.5, 10, by=.1) 
MSE <- matrix(0, length(THETA), 2) 
for(j in 1:length(THETA)) 
{
   D <- matrix(rnorm(k*n, mean=THETA[j], sd=THETA[j]), k, n) 
   ThetaHat_1 <- apply(D, 1, mean)    
   ThetaHat_2 <- apply(D, 1, sd) 
   MSE[j,1] <- mean( (ThetaHat_1 - THETA[j])^2 ) 
   MSE[j,2] <- mean( (ThetaHat_2 - THETA[j])^2 ) 

}   
plot(THETA, MSE[,1], xlab=quote(theta), ylab="MSE", 
main=expression(paste("MSE for each value of ", theta)), 
type="l", col=2, cex.lab=1.3, cex.main=1.5)
lines(THETA, MSE[,2], col=4)
MSE1 <- function(Q) (Q[1]^2)/Q[2] 
MSE2 <- function(Q) 
{
   theta <- Q[1]; n <- Q[2]; 
   G <- gamma(n/2)/gamma( (n-1)/2 )
   bias <- theta * (1 - sqrt(2/(n-1)) * G )
   variance <- (theta^2) * (1 - (2/(n-1)) * G^2 )
   return(bias^2 + variance)
}
THETA <- cbind(matrix( seq(.5, 10, length=100), 100, 1 ), rep(50,100)) 
MSE <- matrix(0, 100, 2) 
MSE[,1] <- apply(THETA, 1, MSE1)
MSE[,2] <- apply(THETA, 1, MSE2) 
plot(THETA[,1], MSE[,1], xlab=quote(theta), ylab="MSE", 
main=expression(paste("MSE for each value of ", theta)), 
type="l", col=2, cex.lab=1.3, cex.main=1.5)
lines(THETA[,1], MSE[,2], col=4)
alpha = .05
k <- 1000 
n <- 10*c(1:5) 
mu_D <- seq(0, 2, by=.1) 
Power <- matrix(0, length(mu_D), 5) 
for(i in 1:5)
{
   for(j in 1:length(mu_D))
   {
      X <- matrix( rnorm(n[i]*k), k, n[i]) 
      Y <- matrix( rnorm(n[i]*k, mean=mu_D[j]), k, n[i]) 
      Xmeans <- apply(X, 1, mean) 
      Ymeans <- apply(Y, 1, mean) 
      T <- sqrt(n[i])*(Xmeans - Ymeans)/sqrt(2) 
      I <- (abs(T) > qnorm(1-(alpha/2))) 
      Power[j,i] <- mean(I)
   }
}      
plot(mu_D, Power[,1], xlab=quote(mu(D)), ylab=expression( 
paste("Power(", mu(D), ")")), col=2, cex.lab=1.3,
cex.main=1.5, main=expression(paste("Power(", mu(D), ") vs.", mu(D))),
type="l" ) 
points(mu_D, Power[,1], col=2)
points(mu_D, Power[,2], col=3)
points(mu_D, Power[,3], col=4)
points(mu_D, Power[,4], col=5)
points(mu_D, Power[,5], col=6)
lines(mu_D, Power[,2], col=3)
lines(mu_D, Power[,3], col=4)
lines(mu_D, Power[,4], col=5)
lines(mu_D, Power[,5], col=6)
legend(1.5, .3, c("n = 10", "n = 20", "n = 30", "n = 40", "n = 50"), pch=(1), 
col=c(2:6), lty=1)
abline(h=alpha)
alpha <- .05
n <- 10*c(1:5) 
mu_D <- seq(0, 2, by=.1) 
Power <- matrix(0, length(mu_D), 5) 
for(i in 1:5)
{
   for(j in 1:length(mu_D))
   {
      Power[j,i] <- 1 - ( pnorm( qnorm(1-alpha/2) - sqrt(n[i])*mu_D[j]/sqrt(2) ) - 
      pnorm( qnorm(alpha/2) - sqrt(n[i])*mu_D[j]/sqrt(2) ) ) 
   }
}
plot(mu_D, Power[,1], xlab=quote(mu(D)), ylab=expression( 
paste("Power(", mu(D), ")")), col=2, cex.lab=1.3,
cex.main=1.5, main=expression(paste("Power(", mu(D), ") vs.", mu(D))),
type="l" ) 
points(mu_D, Power[,1], col=2)
points(mu_D, Power[,2], col=3)
points(mu_D, Power[,3], col=4)
points(mu_D, Power[,4], col=5)
points(mu_D, Power[,5], col=6)
lines(mu_D, Power[,2], col=3)
lines(mu_D, Power[,3], col=4)
lines(mu_D, Power[,4], col=5)
lines(mu_D, Power[,5], col=6)
legend(1.5, .3, c("n = 10", "n = 20", "n = 30", "n = 40", "n = 50"), 
	pch=(1), col=c(2:6), lty=1)
abline(h=.05)
