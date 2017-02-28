## This was the final accept-reject example from Ch 2
f=function(theta){(1+theta^2/2)^-1.5*exp(-(1-theta)^2/2)}
curve(f(x), from=-3, to =4)
g=function(theta){dt(theta-.5, 2)}
ratiofunc=function(theta){f(theta)/g(theta)}
curve(f(x)/g(x), from=-3, to =3)
optimize(ratiofunc,c(-5,5), maximum=T)
M=optimize(ratiofunc,c(0,1), maximum=T)$objective
curve(f(x), from=-3, to =4)
lines(seq(-3,4,length=200),M*g(seq(-3,4,length=200)),col=2)

Nsim=10^4
Y=rt(M*Nsim,2)+.5
U=runif(M*Nsim, max=M)
accept=U<ratiofunc(Y)
posttheta=Y[accept]
mean(posttheta)
hist(posttheta, freq=F)
plot(Y[1:1000], U[1:1000], col=accept[1:1000]+1, pch=16, xlim=c(-4,4),ylim=c(0,2.5))
plot(density(posttheta))
lines(density(Y))

## Using the integrate command for numerical integration
be=function(t,a,b){t^(a-1)*(1-t)^(b-1)}
betaint=function(a,b){integrate(be,0,1,a=a, b=b)$val}
betaout=mapply(betaint, a=rep(1,20), b=seq(.5,5,length=20))
truebeta=mapply(beta, a=rep(1,20), b=seq(.5,5,length=20))
plot(betaout,truebeta)
betaout-truebeta

## Integrate computing a marginal density poorly
## x~N(theta,8^2), theta|x~N(x, 6^2) true marginal of two xs is product of dnorm(x,0,10)

sigmatheta=6
sigmax=8
sigmamarg=sqrt(sigmatheta^2+sigmax^2)
theta=rnorm(1,0,sigmatheta)
xobs=rnorm(2,theta,sigmax)
truemargx=prod(dnorm(xobs,0,sigmamarg))
jointlik=function(the,xobs){prod(dnorm(xobs,the,sigmax))*dnorm(the,0,sigmatheta)}
nummargx=integrate(jointlik, -Inf,Inf,xobs=xobs)
print(c(xobs,truemargx))
print(nummargx)

#Clasical Monte Carlo integration - x^x and some other interesting toy functions
h=function(x){x^x}  #function to be integrated
curve(h,xlab="x", ylab="x^x")
#h=function(x){.75-(.5-x)^2+.5*(1-x)*sin(20*pi*x)}
#h=function(x){1.5-(.5-x)^2+6*(.5-x)^2*sin(20*pi*x), ylim=c(0,3)}
integrate(h,0,1)  #numerical integration estimate
x=h(runif(10^4))  # evaluate h at draws from uniform
estint=cumsum(x)/(1:10^4) #compute the sequence of averages
esterr=sqrt(cumsum((x-estint)^2))/(1:10^4)
plot(estint, xlab="Mean and error range", type="l", lwd=2, ylim=c(mean(x)-20*esterr[10^4],mean(x)+20*esterr[10^4]), ylab="")
lines(estint+1.96*esterr, col="grey", lwd=2)
lines(estint-1.96*esterr, col="grey", lwd=2)


## Binomial tail probability
p=1-pbinom(9,16,.5)
NSim=10^4
x=rbinom(NSim, 16, .5)
h=(x>9)
estint=cumsum(h)/(1:NSim) #compute the sequence of averages
esterr=sqrt(cumsum((h-estint)^2))/(1:NSim)
plot(estint, xlab="Mean and error range", type="l", lwd=2, ylim=c(mean(h)-20*esterr[NSim],mean(h)+20*esterr[NSim]), ylab="")
lines(estint+1.96*esterr, col="grey", lwd=2)
lines(estint-1.96*esterr, col="grey", lwd=2)
lines(p+1.96*sqrt(p*(1-p)/(1:NSim)),col="blue")
lines(p-1.96*sqrt(p*(1-p)/(1:NSim)),col="blue")
lines((1:NSim),rep(p,NSim))

##Importance Sampling - normal from double exponential
curve(.5*dexp(abs(x)),from=-4, to =4)
lines(seq(-4,4,by=.01), dnorm(seq(-4,4,by=.01),0,2),col=2)
Nsim=10^4
#library(VGAM)
#x=rlaplace(Nsim)
x=rexp(Nsim)*(2*rbinom(Nsim,1,.5)-1)
plot(density(x))
w=dnorm(x,0,2)/(.5*dexp(abs(x)))
lines(density(x, weights=w/sum(w)),col=2)

###Importance Sampling - tail probability of t
df=8
v=10
par(mfrow=c(1,1))
truevalue=1-pt(v,df)
Nsim=10^4
x=rexp(Nsim)+v
w=dt(x,df)/dexp(x-v)
plot(cumsum(w)/1:Nsim, type="l", ylim=c(0, truevalue*1.1))
abline(a=truevalue, b=0, col=2)
y=rt(Nsim,df)
lines(cumsum(y>v)/1:Nsim, type="l", col=3)

## Importance sampling for tail area of normal distribution
par(mfrow=c(1,1))
v=5
truevalue=1-pnorm(v)
Nsim=10^5
x=rexp(Nsim)+v
w=dnorm(x)/dexp(x-v)
plot(cumsum(w)/1:Nsim, type="l")
abline(a=truevalue, b=0, col=2)
y=rnorm(Nsim)
lines(cumsum(y>v)/1:Nsim, type="l", col=3)

## Resampling
Nsim=10^4
x=rexp(Nsim)*(2*rbinom(Nsim,1,.5)-1)
w=dnorm(x,0,2)/(.5*dexp(abs(x)))
y=sample(x,Nsim, prob=w, replace=T)
plot(density(x),lwd=3)
lines(density(x, weights=w/sum(w)),col=2,lwd=3)
lines(density(y),col=3,lwd=3)

### Self-normalized weighting for posterior mean
f=function(x){(1+x^2/2)^(-1.5)*exp(-.5*(1-x)^2)}
curve(f(x),-2,4)
Nsim=10^4
x=rt(Nsim,2)+.4
w=f(x)/dt(x-.4,2)
w=w/sum(w)
plot(density(x),xlim=c(-4,4),ylim=c(0,.75))
lines(density(x,weights=w),col=2)
myseq=seq(-10,10,by=.1)
lines(myseq,f(myseq),col=3)
postmean=sum(x*w)
postmean

## Importance Sampling for normal/t Bayes model
f=function(x){(1+x^2/2)^(-1.5)*exp(-.5*(1-x)^2)}
curve(f(x),-2,4)
Nsim=10^4
x=rt(Nsim,2)+.4
w=f(x)/dt(x-.4,2)
w=w/sum(w)
plot(density(x),xlim=c(-4,4),ylim=c(0,.75))
lines(density(x,weights=w),col=2)
myseq=seq(-10,10,by=.1)
lines(myseq,f(myseq),col=3)
postmean=sum(x*w)
postmean

## Bad idea to use normal as instrumental for cauchy
curve(dcauchy(x)/dnorm(x),-5,5)
Nsim=10^6
x=rnorm(Nsim)
w=dcauchy(x)/dnorm(x)
boxplot(w/sum(w))
plot(cumsum(w*x)/cumsum(w),type="l")
abline(a=0,b=0,col=2)