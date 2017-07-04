#uniform draw
Nsim=10000
x=runif(Nsim)

par(mfrow=c(1,3)) 
hist(x)
plot(x[-Nsim],x[-1])
acf(x)


#To implement the inverse transform method, we'll generate from $\mathcal{U}(0,1)$ and take the square root to get values from $Beta(2,1)$
Nsim=10^4
U=runif(Nsim)
X=sqrt(U)
Y=rbeta(Nsim,2,1)
par(mfrow=c(1,2)) 
hist(X, freq=F, main="Beta from Uniform")
hist(Y, freq=F, main="Beta from R")


#Cauchy Inverse Transform
Finv=function(y,mu,sig){mu+sig*tan(pi*(y-.5))}
mu=0
sig=1
Nsim=10^4
U=runif(Nsim)
x=Finv(U,mu,sig)
z=rcauchy(Nsim)
par(mfrow=c(1,2))
hist(x, breaks=c(min(x), seq(-20,20,by=.5), max(x)), freq=F,plot=T, ylim=c(0,.35), xlim=c(-20,20))
hist(z, breaks=c(min(z), seq(-20,20,by=.5), max(z)), freq=F,plot=T, ylim=c(0,.35), xlim=c(-20,20))

#If you already have a standard normal random variable, these values can be squared to give a sample from a chi-squared distribution with 1 degree of freedom (or squared and summed for a chi-square with more degrees of freedom)

Nsim=10^6
system.time(rnorm(Nsim)^2)
system.time(rchisq(Nsim, 1))
X=rnorm(Nsim)^2
Z=rchisq(Nsim, 1)

par(mfrow=c(1,2))
hist(X,  freq=F,plot=T)
hist(Z, freq=F,plot=T)

#The following R code lets us generate from a Binomial(n=12, p=.4).
n=12
p=.4
Nsim=10^4
cump=pbinom(0:n, n, p)
U=runif(Nsim)
x=findInterval(U, cump)
hist(x, freq=F, breaks=c(0:n))

#negative binom 
  cump=pnbinom(c(1:350),50,.2)
round(cump,6)
cumpsubset=cump[102:338]
Nsim=10^4
U=runif(Nsim)
X=findInterval(U, cumpsubset)+102
table(X)/Nsim

## Accept Reject
Nsim=2*10^4
Y=runif(Nsim)
U=runif(Nsim)
accept=(U<1-Y^2)
X=Y[accept]
sum(accept)/Nsim
par(mfrow=c(2,1))
hist(X, freq=F)
plot(Y[1:1000], 1.5*U[1:1000], col=accept[1:1000]+1, pch=16)