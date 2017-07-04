##Exercise 2.6.  Compare execution times of Accept-Reject algorithms.  Variables to compare include 1) should we draw from uniform(0,1) and mulitply by M or draw from uniform(0,M)? and 2) how should we draw until we get the desired number of simulations?

set.seed(907104)
## Set-up for beta accept-reject problem
a=2.7
b=6.3
f=function(x,a,b){dbeta(x,a,b)}
M=optimize(f,c(0,1), a=a, b=b,maximum=T)$objective
Nsim=2500

## Focus on the command of interest, generating from uniform to detect time differences
## Upped number of simulations so time was nonzero
system.time(runif(10^6)*M)
system.time(runif(10^6, max=M))

##Generate Nsim * 1/Acceptance Rate.  Not guaranteed to get Nsim acceptances

u=runif(Nsim*M, max=M)
y=runif(Nsim*M)
mysample=y[u<dbeta(y,a,b)]

accept=u<dbeta(y,a,b)  #Indicator of whether a value was accepted
acceptrate=mean(accept)   ## The acceptance rate.  Should be close to 1/M
length(mysample)         #Number of acceptances
plot(density(mysample), lwd=3, col=3)  #Density of sample
myseq=seq(0,1,by=.01)    #x values for plotting true density
lines(myseq, dbeta(myseq,a,b), lwd=3) #True density
points(y,u,col=accept+1)  #Plot points color coded by whether they were accepted




#Method 2 with while as in book  The time of this method varies considerably depending on whether you have to enter the loop a second time
Bookwhile=function(Nsim, M, a,b){
mysample=NULL
while(length(mysample)<Nsim){
  y=runif(Nsim*M)
  mysample=c(mysample,y[runif(Nsim*M, max=M)<dbeta(y,a,b)])
}
return(mysample[1:Nsim])
}
Bookwhiletimes=replicate(100,system.time(Bookwhile(Nsim,M,a,b)))
mean(Bookwhiletimes["elapsed",])

#While using higher Nprop so second loop less likely to be required
Multiplierwhile=function(Nsim,M,a,b,mult=1.1){
mysample=NULL
while(length(mysample)<Nsim){
  y=runif(Nsim*M*mult)
  mysample=c(mysample,y[runif(Nsim*M*mult, max=M)<dbeta(y,a,b)])
}
return(mysample[1:Nsim])
}

Multwhiletimes=replicate(100,system.time(Multiplierwhile(Nsim,M,a,b)))
mean(Multwhiletimes["elapsed",])

#while- number of draws depending on how many still needed.
Adjustwhile=function(Nsim,M,a,b,mult=1.1){
mysample=NULL
while(length(mysample)<Nsim){
  y=runif((Nsim-length(mysample))*M*mult)
  mysample=c(mysample,y[runif((Nsim-length(mysample))*M*mult, max=M)<dbeta(y,a,b)])
}
return(mysample[1:Nsim])
}
Adjustwhiletimes=replicate(100,system.time(Adjustwhile(Nsim,M,a,b)))
mean(Adjustwhiletimes["elapsed",])
