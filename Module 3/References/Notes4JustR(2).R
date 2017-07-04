matrixerror
  set.seed(921124)
Nsim=10^3
Nparallel=500
h=function(x){x^x}  #function to be integrated
curve(h,xlab="x", ylab="x^x")
x=matrix(h(runif(Nsim*Nparallel)),ncol=Nparallel)  # evaluate h at draws from uniform
estint=apply(x,2,cumsum)/(1:Nsim) #compute the sequence of averages
esterr=sqrt(cumsum((x[,1]-estint[Nsim,1])^2))/(1:Nsim)
plot(estint[,1], xlab="Mean and error range", ylim=c(.7,.85),type="l", lwd=2, ylab="")
y=apply(estint,1,quantile,c(.025,.975))        
polygon(c(1:Nsim, Nsim:1),c(y[1,],rev(y[2,])),col="grey")
lines(estint[,1]+1.96*esterr, col="blue", lwd=2)
lines(estint[,1]-1.96*esterr, col="blue", lwd=2)

##bootstraperror
boot=matrix(sample(x[1,],Nparallel*Nsim,rep=T),nrow=Nsim, ncol=Nparallel)
bootest=apply(boot,2,cumsum)/(1:Nsim)
bootupper=apply(bootest,1,quantile, .975)
bootlower=apply(bootest,1,quantile, .025)
lines(bootupper, lty=3)
lines(bootlower, lty=3)

## Importance Sampling Error
f=function(x){(1+x^2/2)^(-1.5)*exp(-.5*(1-x)^2)}
#curve(f(x),-2,4)
Nsim=10^3
x=rt(Nsim,2)+.4
w=f(x)/dt(x-.4,2)
w=w/sum(w)
#plot(density(x),xlim=c(-4,4),ylim=c(0,.75))
#lines(density(x,weights=w),col=2)
#myseq=seq(-10,10,by=.1)
#lines(myseq,f(myseq),col=3)
postmean=cumsum(x*w)/cumsum(w)
plot(postmean, type="l", ylim=c(0,1))
est=postmean[Nsim]
err=sum(w*(x-est)^2)/Nsim*sum(w)*(1+Nsim^2*var(w)/(sum(w))^2) ##MC estimate of error for importance sampling estimate
err

## Effective Sample Size
set.seed(9241141)
Nsim=10^4
isess=function(df, shift){
  x=rt(Nsim,df)+shift
  w=f(x)/dt(x-shift,df)
  w=w/sum(w)
  return(1/sum(w^2))
}

isess(2,0)
isess(2,.4)

## Perplexity
set.seed(9241141)
Nsim=10^4
isperp=function(df, shift){
  x=rt(Nsim,df)+shift
  w=f(x)/dt(x-shift,df)
  wn=pmax(w/sum(w), rep(10^-10,Nsim)) ## Avoid zeros by assigning a minimum weight
  wn=wn/sum(wn)
  return(exp(-sum(wn*log(wn)))/Nsim)
}

isperp(2,0)
isperp(2,.4)

##Rao-Blackwell
Nsim=10^3
a=8
b=.5
trueval=b/(a-1)
y=matrix(1/rgamma(100*Nsim,a, rate=b),ncol=100)  ## generate variances with inversegamma(shape=a, scale=b)
x=matrix(rnorm(100*Nsim,sd=sqrt(y)),ncol=100) ## generate x values using sqrt(y) as standard deviation

curve(sqrt(a/b)*dt(x*sqrt(a/b),df=2*a),-4,4,ylim=c(0,1.6))
lines(density(x),col=2)

matplot(apply(x^2,2,cumsum)/(1:Nsim),type="l",col="grey80",
        lty=1,ylim=trueval*c(.2,2), xlab="",ylab="") ##Find sequences of estimates of x^2 and plot them in light gray
matplot(apply(y,2,cumsum)/(1:Nsim),type="l",col="grey40",
        lty=1,add=T,xlab="",ylab="") ## The RB estimate of x^2 is y.  Plot these sequences in dark gray
abline(h=trueval,col="gold",lty=2,lwd=2) ## The true expected value is b/(a-1)

tradest=apply(x^2,2,cumsum)/(1:Nsim)
rbest=apply(y,2,cumsum)/(1:Nsim)
errfun=function(x){(x-trueval)^2}
errfun2=function(x){sqrt(cumsum((x-trueval)^2))/(1:Nsim)}
traderr=apply(apply(x^2,2,errfun2),1,mean)
rberr=apply(apply(y,2,errfun2),1,mean)
#rberr=apply(sqrt(errfun(rbest))/1:Nsim,1,mean)
plot(traderr/rberr)
traderr[1:20]/rberr[1:20]


##HalfWidth
  set.seed(101121041)
Nsim=2000
x=rnorm(Nsim)
estint=cumsum(x)/1:Nsim
estvar=cumsum((x-estint)^2)/1:Nsim
halfalpha=1.96*sqrt(estvar/1:Nsim)
cbind(1400:1600, halfalpha[1400:1600])
sum(halfalpha>.05)  

#Antithetic
  set.seed(930315)
Nsim=10^3
u=runif(Nsim)
ua=as.vector(rbind(u[1:(Nsim/2)],1-u[1:(Nsim/2)])) #Antithetic set-up
finv=function(y,mu,b){mu-b*log(1/y-1)} # Inverse cdf for logistic
mu<-0 #mean
b<-1  #scale
x<-finv(u,mu,b)
xa=finv(ua,mu,b)
## Functions for cumsum and variance of one sequence
estintfun=function(x){cumsum(x)/1:length(x)}
esterrfun=function(x)
{
  estint=estintfun(x)
  return(sqrt(cumsum((x-estint)^2))/(1:length(x)))
}

##Get sequence of estimates for regular and antithetic uniforms
estint=estintfun(x)
esterr=esterrfun(x)
estinta=estintfun(xa)
esterra=esterrfun(xa)

## Plot
plot(estint,type="l", col="blue", lwd=2, ylim=c(mu-1.5, mu+1.5))
abline(h=mu, col="black", lwd=3)  ## Truth is known
lines(estint+1.96*esterr, col="lightblue", lwd=1) ## Confidence bands based on one sequence
lines(estint-1.96*esterr, col="lightblue", lwd=1)
lines(estinta, col=2, lwd=2)
lines(estinta+1.96*esterra, col="pink", lwd=1) ## Confidence bands based on one sequence
lines(estinta-1.96*esterra, col="pink", lwd=1)