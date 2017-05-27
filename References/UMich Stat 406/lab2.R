 pnorm(1) - pnorm(0)
X <- rnorm(10000)
mean(X <= 1) - mean(X <= 0) 
pexp(1) - pexp(0)
X <- rexp(10000)
mean(X <= 1) - mean(X <= 0) 
punif(1,a=0,b=3) - punif(0,a=0,b=3)
X <- runif(10000,a=0,b=3)
mean(X <= 1) - mean(X <= 0) 
pchisq(1,df=2) - pchisq(0,df=2)
X <- rchisq(10000,df=2)
mean(X <= 1) - mean(X <= 0) 
setwd("users/jasoneg/Desktop")
Y <- read.table('testdata.txt', header=T, row.names=1) 
Y$status 
A <- read.table('flow-occ-table.txt', header=T, sep=",")
V <- scan("testdata2.txt", sep=",") 
X <- list() 
N <- sum( abs(V) > 3 ) 
M <- 0
for(n in 1:N) 
{

   w <- min( which( abs(V) > 3) ) 
   X[[n]] <- V[1:w] 
   V <- V[-c(1:w)]
   if( w > M ) 
   {
      M = w
   }
}
D <- matrix(NA, N, M) 
for(i in 1:N) D[i,1:length(X[[i]])] <- X[[i]]
v <- rnorm(30) 
X <- matrix(v, 10, 3) 
X <- as.data.frame(X) 
X <- cbind(X, c( rep("treatment", 5), rep("control", 5)) ) 
colnames(X) = c("y", "x1", "x2", "status")
write.table(X, file="testdata.txt", row.names=T, col.names=T) 
x1 <- matrix(0, 20, 20) 
x2 <- c("Treatment", "Control") 
x3 <- c(1:10) 
x4 <- matrix(2, 13, 4) 
save("x1","x2","x3","x4",file="variables.Rdata")
load("variables.Rdata")
Traffic <- read.table('flow-occ-table.txt', header=T, sep=",")
colnames(Traffic)
colnames(Traffic) = c("O1", "F1", "O2", "F2", "O3", "F3")
colnames(Traffic)
newTraffic <- matrix(0, nrow(Traffic), 2) 
for(i in 1:nrow(Traffic)) 
{

   traffic.i <- Traffic[i,]
   maxflow <- max( traffic.i[c(2,4,6)] )  
   maxocc <- as.numeric( traffic.i[c(1,3,5[which.max(traffic.i[c(2,4,6)])] )
   newTraffic[i,] <- c(maxflow, maxocc) 

}
newTraffic <- as.data.frame(newTraffic)
colnames(newTraffic) = c("maxflow", "maxocc")
write.table(newTraffic, file="flow-occ-table_clean.txt", sep="\t", header=T, 
   row.names=1)
Y <- read.table('flow-occ-table.txt', sep=",", header=T)
Y <- Y[1:30,]
v <- seq(0, 1, length=10)
plot(v)  
v <- seq(0, 1, length=10)
u <- seq(3, 4, length=10)
plot(u,v) 
plot(Y$Occ1, col=2, ylim=c(0, .04))
points(Y$Occ2, col=3)
points(Y$Occ3, col=4)
plot(Y$Occ1, col=2, ylim=c(0, .04))
points(Y$Occ2, col=3)
points(Y$Occ3, col=4)
plot(Y$Occ1, col=2, ylim=c(0, .04), pch=3)
points(Y$Occ2, col=3, pch=0) 
points(Y$Occ3, col=4, pch=16)
plot(Y$Occ1, col=2, ylim=c(0, .04), pch=3, cex=(2/3))
points(Y$Occ2, col=3, pch=0, cex=(2/3))
points(Y$Occ3, col=4, pch=16, cex=(2/3))
lines(Y$Occ1, col=2)
lines(Y$Occ2, col=3) 
lines(Y$Occ3, col=4) 
plot(Y$Occ1, col=2, ylim=c(0, .04), pch=3, cex=(2/3), xlab="x-axis label", 
   ylab="y-axis label", main="Example plot")
points(Y$Occ2, col=3, pch=0, cex=(2/3))
points(Y$Occ3, col=4, pch=16, cex=(2/3))
lines(Y$Occ1, col=2)
lines(Y$Occ2, col=3) 
lines(Y$Occ3, col=4)
plot(Y$Occ1, col=2, ylim=c(0, .04), pch=3, cex=(2/3), xlab="x-axis label", 
   ylab="y-axis label", main="Example plot", cex.lab=1.25, cex.main=1.5)
points(Y$Occ2, col=3, pch=0, cex=(2/3))
points(Y$Occ3, col=4, pch=16, cex=(2/3))
lines(Y$Occ1, col=2)
lines(Y$Occ2, col=3) 
lines(Y$Occ3, col=4)
plot(Y$Occ1, col=2, ylim=c(0, .04), pch=3, cex=(2/3), xlab="x-axis label", 
   ylab="y-axis label", main="Example plot", , cex.lab=1.25, cex.main=1.5, axes=FALSE)
points(Y$Occ2, col=3, pch=0, cex=(2/3))
points(Y$Occ3, col=4, pch=16, cex=(2/3))
lines(Y$Occ1, col=2)
lines(Y$Occ2, col=3) 
lines(Y$Occ3, col=4)
axis(1, seq(0, 30, length=7), c("a", "b", "c", "d", "e", "f", "g"))
axis(2, seq(0, .04, length=5), c("very low", "low", "medium", "high", "very high"))
box()
plot(Y$Occ1, col=2, ylim=c(0, .04), pch=3, cex=(2/3), xlab="x-axis label", 
   ylab="y-axis label", main="Example plot", , cex.lab=1.25, cex.main=1.5, axes=FALSE)
points(Y$Occ2, col=3, pch=0, cex=(2/3))
points(Y$Occ3, col=4, pch=16, cex=(2/3))
lines(Y$Occ1, col=2)
lines(Y$Occ2, col=3) 
lines(Y$Occ3, col=4)
axis(1, seq(0, 30, length=7), c("a", "b", "c", "d", "e", "f", "g"))
axis(2, seq(0, .04, length=5), c("very low", "low", "medium", "high", "very high"))
box()
legend(25, .035, c("Occ1", "Occ2", "Occ3"), pch=c(3,0,16), col=c(2:4), lty=1)
text(15, .035, "text", cex=5, col=7)
plot.env <- par(bg="green", mfrow=c(2,2), cex.main=1.125, cex.lab= 1.25, 
   col.axis="red", col.lab="blue", col.main="black")
plot( rnorm(100), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(0,1) variables")
plot( rnorm(100,mean=.5), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(.5,1) variables")
plot( rnorm(100,mean=1), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(1,1) variables")
plot( rnorm(100,mean=1.5), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(1.5,1) variables")
pdf("plot.pdf")  
plot.env <- par(bg="green", mfrow=c(2,2), cex.main=1.125, cex.lab= 1.25, 
   col.axis="red", col.lab="blue", col.main="black")
plot( rnorm(100), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(0,1) variables")
plot( rnorm(100,mean=.5), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(.5,1) variables")
plot( rnorm(100,mean=1), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(1,1) variables")
plot( rnorm(100,mean=1.5), xlab="Index", ylab="y-axis", 
   main="Scatterplot of 100 N(1.5,1) variables") 
dev.off() 
