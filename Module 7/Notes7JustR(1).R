## A very basic two-stage Gibbs sampler (for which a chain is really unnecessary)
mu=10
Nsim=10^4
theta=rnorm(1,mu,1)
x=numeric()
for(i in 1:Nsim)
{
  x[i]=rnorm(1,theta[i],1)
  theta=c(theta,rnorm(1,mu,1))  ##more efficient to do this outside loop
}
curve(dnorm(x,mu,sqrt(2)),xlim=c(6,14))
lines(density(x),col=2)


### Bivariate Normal
library(MASS)
Nsim=10^4
rho=.9
#x=mvrnorm(Nsim, c(0,0), matrix(nrow=2,ncol=2,data=c(1,rho,rho,1)))
x1=as.vector(rnorm(1))
x2=numeric(length=Nsim)
csd=sqrt(1-rho^2)  #conditional standard deviation
for(i in 1:Nsim)
{
  x2[i]=rnorm(1,rho*x1[i],csd)
  x1[i+1]=rnorm(1,rho*x2[i],csd)
}
d=kde2d(x1[1:Nsim],x2)
persp(d)
image(d)

###Normal with unknown mean and unknown variance
set.seed(701)
n=30
truetheta=50
truesigma2=5
Nsim=10^4
burnin=100
x=rnorm(n,truetheta,sqrt(truesigma2))
xbar=mean(x)
xvar=var(x)
theta=numeric(length=Nsim)
theta[1]=xbar+rnorm(1,0,.5)
tausq=numeric(length=Nsim) ##Precision
for(i in 1:Nsim)
{
  tausq[i]=rgamma(1,n/2,((n-1)*xvar+n*(xbar-theta[i])^2)/2 )
  theta[i+1]=rnorm(1,xbar,sqrt(1/(n*tausq[i])))
}
plot(density(theta[burnin:Nsim]),col=3)

plot(density(1/tausq[burnin:Nsim]))
plot(density(tausq[burnin:Nsim]))
myseq=seq(.05,10, by=.01)
lines(myseq,(n-1)*xvar*dchisq(myseq*xvar*(n-1),n-1),col=2)


### Gibbs Sampler for Changepoint in Poisson
coal.mining.disasters<-c(4,5,4,1,0,4,3,4,0,6,3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,
                         4,2,5,2,2,3,4,2,1,3,2,2,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,0,1,1,1,0,1,0,1,0,0,
                         0,2,1,0,0,0,1,1,0,2,3,3,1,1,2,1,1,1,1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1)

kprob<-function(k,y,lambda,phi,a, b, c, d, n)
{
  if(k<n)
  {
    return(lambda^(a-1+sum(y[1:k]))*phi^(c-1+sum(y[(k+1):n]))*
      exp(-(k*lambda+(n-k)*phi+b*lambda+d*phi)))
  }
  else{return(lambda^(a-1+sum(y[1:k]))*phi^(c-1)*exp(-(k*lambda+b*lambda+d*phi)))}
}

gibbs.poisson.gamma<-function(theta.matrix, y, reps)
{
  a<-4;b<-1; c<-1; d<-2
  n<-length(y)
  kseq<-c(1:n)
  for(i in 2:(reps+1))
  {
    lambda<-rgamma(1,a+sum(y[1:theta.matrix[(i-1),3]]),(b+theta.matrix[(i-1),3]))
    phi<-rgamma(1,c+sum(y[(theta.matrix[(i-1),3]+1):n]),(d+n-theta.matrix[(i-1),3]))
    k<-sample(kseq,1, prob=lapply(kseq,FUN=kprob,y,lambda, phi, a, b, c, d, n))
    theta.matrix<-rbind(theta.matrix, c(lambda, phi, k))
  }
  return(theta.matrix )
}

Nsim=10^4
start<-matrix(c(1,1,40),1,3)
#coal.sample<-gibbs.poisson.gamma(start,coal.mining.disasters, Nsim)
coal.sample<-gibbs.poisson.gamma(start, coal.mining.disasters, 1000)



plot.walk.lambda<-function(walk.mat)
{
  walk.mat<-walk.mat[,-2]
  plot(walk.mat[1,1], walk.mat[1,2], xlim=c(0,5), ylim=c(15,60), xlab="lambda", ylab="k")
  for(i in 1:(nrow(walk.mat)-1)) 
  {
    segments(walk.mat[i,1], walk.mat[i,2],walk.mat[(i+1),1],walk.mat[i,2])
    segments(walk.mat[(i+1),1], walk.mat[i,2], walk.mat[(i+1),1], walk.mat[(i+1),2])
  }
  
}


plot.walk.phi<-function(walk.mat)
{
  walk.mat<-walk.mat[,-1]
  plot(walk.mat[1,1], walk.mat[1,2], xlim=c(0,2), ylim=c(15,60), xlab="phi", ylab="k")
  for(i in 1:(nrow(walk.mat)-1)) 
  {
    segments(walk.mat[i,1], walk.mat[i,2],walk.mat[(i+1),1],walk.mat[i,2])
    segments(walk.mat[(i+1),1], walk.mat[i,2], walk.mat[(i+1),1], walk.mat[(i+1),2])
  }
  
}


plot.walk.lambda(coal.sample[1:100,])
plot.walk.lambda(coal.sample[101:201,])

plot.walk.phi(coal.sample[1:100,])

plot(density(coal.sample[,2]),xlim=c(0,5))
lines(density(coal.sample[,1]),col=2)

plot(hist(coal.sample[,3], breaks=c(30:50)))

##Truncated Exponential
Nsim=10^4
set.seed(846)
truetheta=.1
a=9 ##censoring time
n=20
xu=rexp(n,truetheta)
## Determine observed vs. missing data
y=replace(xu,xu>a,a)
y=sort(y)
m=sum(xu<a)
x=y[1:m]
theta=numeric()
xbar=mean(x)
theta[1]=1/xbar
z=matrix(nrow=Nsim, ncol=n-m)
for(i in 1:Nsim)
{
  z[i,]=rep(a,n-m)+rexp(n-m,theta[i])
  theta[i+1]=rgamma(1,n+1,m*xbar+sum(z[i,]))
}
plot(density(theta))
mean(theta)
thetaMLE=1/(xbar+(n-m)*a/m)
quantile(theta,c(.025, .975))

## Comparison of two random effects model parameterizations.  Install the library mcsm first
library(mcsm)
data(Energy)
Energy
randomeff
reparareff
reffout1=randomeff(5.5*10^3,10,30)
reffout2=reparareff(5.5*10^3,10,30)

### Rao-Blackwell Precision 
set.seed(701)
n=30
truetheta=50
truesigma2=5
Nsim=10^4
x=rnorm(n,truetheta,sqrt(truesigma2))
xbar=mean(x)
xvar=var(x)
theta=numeric(length=Nsim)
theta[1]=xbar+rnorm(1,0,.5)
tausq=numeric(length=Nsim) ##Precision
etausq=numeric(length=Nsim)
for(i in 1:Nsim)
{
  tausq[i]=rgamma(1,n/2,((n-1)*xvar+n*(xbar-theta[i])^2)/2 )
  etausq[i]=n/((n-1)*xvar+n*(xbar-theta[i])^2)
  theta[i+1]=rnorm(1,xbar,sqrt(1/(n*tausq[i])))
}

plot(density(etausq),col=3, xlim=c(0,.4))
lines(density(tausq))
mean(etausq)
mean(tausq)
var(etausq)
var(tausq)
var(tausq)/var(etausq)
plot(cumsum(tausq)/(1:Nsim),type="l")
lines(cumsum(etausq)/(1:Nsim),col=3)


## Blocked Gibbs Sampler
library(MASS)
unblock=function(Nsim,a){
  x=as.vector(rnorm(1))
  y=as.vector(rnorm(1))
  z=numeric(Nsim)
  sdyz=sqrt(-3/(4*(a-1)*(a+1)))
  sdx=sqrt((a-1)*(a+1)/(a^2-a/2-.5))
  for(i in 1:Nsim)
  {
    z[i]=rnorm(1,-(2/3)*((a-1)*x[i]-(.5-2*a)*y[i]),sdyz)
    x[i+1]=rnorm(1,(y[i]+z[i])/(2*(a+1)),sdx)
    y[i+1]=rnorm(1,-(2/3)*((a-1)*x[i+1]-(.5-2*a)*z[i]),sdyz)
  }
  return(cbind(x[1:Nsim],y[1:Nsim],z))
}
Nsim=10^4
a=.75
xyzout=unblock(Nsim,a)
# d=kde2d(xyzout[,1],xyzout[,2])
# persp(d)
# image(d)

blockmn=function(Nsim,a){
  x=as.vector(rnorm(1))
  yz=matrix(nrow=Nsim, ncol=2)
  varmat=solve(matrix(nrow=2,ncol=2,data=c(.75, a-.25, a-.25, .75)))
  for(i in 1:Nsim)
  {
    yz[i,]=mvrnorm(1,c(x[i]/2,x[i]/2),varmat)
    x[i+1]=rnorm(1,(yz[i,1]+yz[i,2])/(2*(a+1)),sqrt((a-1)*(a+1)/(a^2-a/2-.5)))
  }
  return(cbind(x[1:Nsim],yz))
}

xyzoutb=blockmn(Nsim,a)
acf(xyzout[,1])
acf(xyzoutb[,1])

# system.time(unblock(Nsim,a=.2))
# system.time(blockmn(Nsim,a=.2))
# 
# system.time(unblock(Nsim,a=.9))
# system.time(blockmn(Nsim,a=.9))