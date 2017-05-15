

#Harry
P=matrix(nrow=3, ncol=3, data=c(0,1/3,1,1/3,0,0,2/3,2/3,0))
I=diag(3)
spi=c(.45, .15, .4)
t(P)%*%spi
#stapi=solve(t(P)-I,matrix(ncol=1,nrow=3, data=rep(0,3)))
Nsim=10^4
x=sample(1:3,1,replace=F,prob=spi)
for(i in 1:Nsim)
  {x=c(x,sample(1:3,1,replace=F, prob=P[x[length(x)],]))}
plot(x[1:25], type="l")
table(x)
#Excercise 6.1
rho=.9
Nsim=10^4
x=rnorm(1)
for(i in 1:Nsim)
{
  x=c(x, rho*x[i]+rnorm(1))
}
plot(density(x))
myseq=seq(-7,7,length=1001)
lines(myseq, dnorm(myseq,0,sqrt(1/(1-rho^2))),col=2)
rho=1
x=rnorm(1)
for(i in 1:Nsim)
{
  x=c(x, rho*x[i]+rnorm(1))
}
plot(density(x))
#Beta
a=4
b=6
Nsim=10^4
X=rep(runif(1), Nsim) ## Initialize chain and set up storage
for(i in 2:Nsim)
  {
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X[i-1],a,b)
  X[i]=X[i-1]+(Y-X[i-1])*(runif(1)<rho)
  }
plot(density(X))
lines(seq(0,1,length=101), dbeta(seq(0,1,length=101),a,b),col=2)
acf(X)

#MHRaschSimple
f<-function(theta,y,b){prod(dbinom(y,1,1/(1+exp(-(theta-b)))))*exp(-theta^2/2)}

MHRasch<-function(y,b,N,theta.start)
{
theta<-c(theta.start)
for(i in 2:(N+1))
{
thetanew<-rnorm(1,theta[i-1],1)
a<-min(1, f(thetanew,y,b)/f(theta[i-1],y,b))
u<-runif(1,0,1)
if(a>u){theta<-c(theta,thetanew)}
else{theta<-c(theta, theta[i-1])}
}
return(theta)
}

theta.start<-c(0)
y<-c(1,1,1,1,0)
b<-c(-1,0,-2,1,2.5)
theta<-MHRasch(y,b,5000,theta.start)

plot(density(theta))
mean(theta)
sd(theta)

#MHRaschFull
f<-function(theta,y,b){prod(dbinom(y,1,1/(1+exp(-(theta-b)))))*exp(-theta^2/2)}

MHRasch<-function(y,b,N,theta.start)
{
myseq<-seq(-3,3, length=6001)
theta<-c(theta.start)
for(i in 2:(N+1))
{
thetanew<-rnorm(1,theta[i-1],1)
thetalast<-theta[length(theta)]
plot(myseq, mapply(f, theta=myseq, MoreArgs=list(y=y,b=b)), type="l",xlab=expression(theta),ylab="", yaxt="n")
points(thetanew, f(thetanew,y,b), pch=16, col=2)
points(thetalast, f(thetalast,y,b), pch=16)
a<-min(1, f(thetanew,y,b)/f(thetalast,y,b))
u<-runif(1,0,1)
print(c(thetanew, thetalast,f(thetanew,y,b), f(thetalast,y,b),a,u))
if(a<1)
{
lines(c(thetanew, thetanew), c(0,f(thetalast,y,b)))
points(thetanew, u*f(thetalast,y,b),pch="u")
}
if(u>a){points(thetanew, f(thetanew,y,b), pch=4, cex=2)}
else{points(thetanew, f(thetanew,y,b), pch=1, cex=2)}
if(a>u){theta<-c(theta,thetanew)}
else{theta<-c(theta, thetalast)}
}
return(theta)
}

theta<-c(0)
y<-c(1,1,1,1,0)
b<-c(-1,0,-2,1,2.5)
theta<-MHRasch(y,b,1,theta)

#MHnormmix
truep=.7
n=100
Nsim=10^4
z=rbinom(n,1,truep)
y=rnorm(n,10-3*z,.5)
lik=function(x){prod(x*dnorm(y,7,.5)+(1-x)*dnorm(y,10,.5))}
p=runif(1)

##Uniform proposal.  No q needed in acceptance prob since q=1
for(i in 2:Nsim)
  {
  prop=runif(1)
  move=runif(1)<(lik(prop)/lik(p[i-1]))
  p=c(p,p[i-1]*(1-move)+prop*move)               
  }
plot(density(p), xlim=c(0,1))
mean(z)

##Beta proposal.  This is an example of a poor proposal since it has little probability near the true posterior mean.
p2=runif(1)
for(i in 2:Nsim)
  {
  prop=rbeta(1,2,10)
  move=runif(1)<(lik(prop)*dbeta(p2[i-1],2,10)/(lik(p2[i-1])*dbeta(prop,2,10)))
  p2=c(p2,p2[i-1]*(1-move)+prop*move)               
  }
lines(density(p2),col=2)



##HastingsNorm

MHNorm=function(Nsim, delta)
{
x=runif(1,-delta,delta)
for(i in 2:Nsim)
  {
  prop=x[i-1]+runif(1,-delta,delta)
  move=runif(1)<dnorm(prop)/dnorm(x[i-1])
  x=c(x,prop*move+x[i-1]*(1-move))
  }
return(x)
}

x1=MHNorm(10^4,.1)
x2=MHNorm(10^5,.1)
x3=MHNorm(10^4, 1)
x4=MHNorm(10^5, 1)
x5=MHNorm(10^4, 10)
x6=MHNorm(10^5, 10)

accept=function(x){mean(x[-c(1)]-x[-c(length(x))]!=0)}
accept(x1)
accept(x2)
accept(x3)
accept(x4)
accept(x5)
accept(x6)

curve(dnorm(x), ylim=c(0, .5), xlim=c(-3.5, 3.5),lwd=3)
lines(density(x1),col=4)
lines(density(x2), col=4,lty=2)
lines(density(x3), col=2)
lines(density(x4), col=2, lty=2)
lines(density(x5), col=3)
lines(density(x6), col=3, lty=2)
#MHnormmixaccept
truep=.7
n=100
Nsim=10^4
z=rbinom(n,1,truep)
y=rnorm(n,10-3*z,.5)
lik=function(x){prod(x*dnorm(y,7,.5)+(1-x)*dnorm(y,10,.5))}
p=runif(1)

##Uniform proposal.  No q needed in acceptance prob since q=1
for(i in 2:Nsim)
  {
  prop=runif(1)
  move=runif(1)<(lik(prop)/lik(p[i-1]))
  p=c(p,p[i-1]*(1-move)+prop*move)               
  }
plot(density(p), xlim=c(0,1))
mean(z)

##Beta proposal.  This is an example of a poor proposal since it has little probability near the true posterior mean.
p2=runif(1)
for(i in 2:Nsim)
  {
  prop=rbeta(1,2,10)
  move=runif(1)<(lik(prop)*dbeta(p2[i-1],2,10)/(lik(p2[i-1])*dbeta(prop,2,10)))
  p2=c(p2,p2[i-1]*(1-move)+prop*move)               
  }
lines(density(p2),col=2)

accept=function(x){mean(x[-c(1)]-x[-c(length(x))]!=0)}
accept(p)
accept(p2)

a=10
b=5
p3=runif(1)
for(i in 2:Nsim)
  {
  prop=rbeta(1,a,b)
  move=runif(1)<(lik(prop)*dbeta(p3[i-1],a,b)/(lik(p2[i-1])*dbeta(prop,a,b)))
  p3=c(p3,p3[i-1]*(1-move)+prop*move)               
  }
lines(density(p3),col=3)

accept(p3)

