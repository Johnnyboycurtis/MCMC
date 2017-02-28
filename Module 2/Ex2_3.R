##Ex 2.3
set.seed(908956)

mymat=matrix(nrow=12, ncol=10000, data=runif(120000,-.5,.5)) ##Uniform values to be summed for CLT approximation
CLTnorm=colSums(mymat)

BoxMmat=matrix(nrow=2,ncol=5000, data=runif(10000)) ##Uniform values to be transformed for Box-Muller
BoxMnorm=sqrt(-2*log(BoxMmat[1,]))*cos(2*pi*BoxMmat[2,])
BoxMnorm=append(BoxMnorm, sqrt(-2*log(BoxMmat[1,]))*sin(2*pi*BoxMmat[2,]))

myseq=seq(-4,4,by=.01)
plot(myseq, dnorm(myseq),type="l", lwd=2)  #True normal density
lines(density(CLTnorm),col=2, lwd=2)       #Smoothed histogram of CLT random variables
lines(density(BoxMnorm),col=3, lwd=2)      #Smoothed histogram of Box-Muller random variables
par(mfrow=c(1,2))
hist(CLTnorm, plot=T, freq=F, xlim=c(-5,-3))  ##Histogram of left tail of CLT
hist(BoxMnorm, plot=T, freq=F, xlim=c(-5,-3))  ##Histogram of left tail of Box-Muller
par(mfrow=c(1,1))

## Look at how the quantiles match up in the tail
normqs=qnorm(seq(.00005,.99995,by=.0001))     ##Normal quantiles    
plot(seq(.00005,.99995,by=.0001), normqs, type="l",xlim=c(0.0,0.01), lwd=2) 
lines(seq(.00005,.99995,by=.0001), sort(BoxMnorm), col=3, lwd=2)
lines(seq(.00005,.99995,by=.0001), sort(CLTnorm), col=2, lwd=2)



## Histograms are not that great.  Instead I have done an estimate of P(z<-2.326348), which should be .01
###  Draw 500 samples with each method and look at how accurately each estimates P(z<-2.326348)
storeCLT=vector()
storeBoxM=vector()
for(i in 1:500){
  mymat=matrix(nrow=12, ncol=10000, data=runif(120000,-.5,.5))
  CLTnorm=apply(mymat,2,sum)
  BoxMmat=matrix(nrow=2,ncol=5000, data=runif(10000))
  BoxMnorm=sqrt(-2*log(BoxMmat[1,]))*cos(2*pi*BoxMmat[2,])
  BoxMnorm=append(BoxMnorm, sqrt(-2*log(BoxMmat[1,]))*sin(2*pi*BoxMmat[2,]))
  storeCLT[i]=sum(CLTnorm<qnorm(.01))/10000
  storeBoxM[i]=sum(BoxMnorm<qnorm(.01))/10000
}
plot(density(storeBoxM), ylim=c(0,500))
lines(density(storeCLT), col=2)
c(mean(storeCLT), sd(storeCLT))
c(mean(storeBoxM), sd(storeBoxM))

