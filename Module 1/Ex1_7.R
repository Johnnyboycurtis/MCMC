## Ex 1.7a
set.seed(829914)
y=c(4.313, 4.513, 5.489, 4.265, 3.641, 5.106, 8.006, 5.087) ##Data set
nBoot=1000      ## How many times to bootstrap
Bmean=replicate(nBoot,mean(sample(y, size=8, replace=T)))
hist(Bmean,breaks=seq(3,8,by=.2), freq=F)   #Frequency is true by default.  Setting it to false gives a relative frequency histogram
lines(density(Bmean),col=2)  #Adds the smoothed curve over the histogram
lines(seq(3,8,by=.1), dnorm(seq(3,8,by=.1), mean(Bmean), sd(Bmean))) #Adds the normal approx over the histogram
quantile(Bmean, .95)  ## Find the 95th percentile of the vector Bmean

## Ex 1.7b
set.seed(618132)
y=c(4.313, 4.513, 5.489, 4.265, 3.641, 5.106, 8.006, 5.087) ##Data set
nBoot1=1000  #number of samples to take from y
nBoot2=1000  #number of times to compute qhat.95 for confidence interval construction
Bq=vector(length=nBoot2) #storage for quantiles
for(i in 1:nBoot2){
  Bmean=replicate(nBoot1,mean(sample(y, size=8, replace=T)))
  Bq[i]=quantile(Bmean, .95)
}
quantile(Bq, c(.025, .975))