### Monte Carlo p-value

## The manufacturer claims the time until failure of a part is gamma(10,1).  10 of the parts get
## installed in your product and your product fails after the 7th failure.
## You cannot observe the failure times of the earlier failures, only when the product fails.
## You observe failure times of 7.46, 4.97, 3.46, 7.17, 9.36, 7.59, 7.19, 8.98, 9.89, 10.90 in 
## your quality control checks.
## a. If the manufacturer's claim is true, what is the probability the average 7th order statistic will
## be as low or lower than observed?
## b. If the manufacturer's claim is true, what the the probability the smallest of 10 7th order statistics
### will be as low or lower than observed?

#failtimes=c(7.46, 4.97, 3.46, 7.17, 9.36, 7.59, 7.19, 8.98, 9.89, 10.90)
failtimes=c(10.46, 8.97, 9.46, 8.17, 9.36, 11.59, 10.19, 8.98, 9.89, 10.90)
Nsim=10000
Nsim=10000
sims=array(dim=c(10,Nsim,10), data=rgamma(100*Nsim,10,1))
simstemp=array(dim=c(10,Nsim,10))
for(i in 1:10){for(j in 1:Nsim){simstemp[i,j,]=sort(sims[i,j,])}}
sims=simstemp[,,7]
observedmean=mean(failtimes)
simmeans=apply(sims,2,mean)
meanpval=mean(simmeans<observedmean)
simmin=apply(sims,2,min)
observedmin=min(failtimes)
minpval=mean(simmin<observedmin)
