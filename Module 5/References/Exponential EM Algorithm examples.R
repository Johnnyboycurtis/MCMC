## Traditional EM for censored exponential
## Generate Data
#set.seed(701)
set.seed(846)
truetheta=10
a=9 ##censoring time
n=20
x=rexp(n,1/truetheta) #R uses rate

## Determine observed vs. missing data
y=replace(x,x>a,a)
y=sort(y)
r=sum(x<a) ## observed data
yc=y[1:r] ## observed data
ycbar=mean(yc) ## mean of the observed data
### EM
thetahat = ycbar ## MLE using observed data
cur = ycbar ## current estimate of theta for EM loop
diff=1 ## set difference variable for while-loop
while(diff>10^-4)
{   
    i = length(thetahat)
    new_estimate = (ycbar*r+(n-r)*(a+thetahat[i]))/n
    thetahat=c(thetahat, new_estimate) 
    diff=abs(thetahat[length(thetahat)]-thetahat[length(thetahat)-1])
}
plot(thetahat, type="l")
thetaMLE=ycbar+(n-r)*a/r
abline(h=thetaMLE, col=2)



mle_f <- function(x){
    "x will be theta"
    out = (3*11+(20-11)*(9+x))/20
    #out = (ycbar*r+(n-r)*(a+x))/n
    return(out)
}

curve(mle_f, from = 0, to = 20)

## Monte Carlo EM for censored exponential
thetahat=cur=ycbar
m=100
diff=1
while(diff>10^-4)
{
  iter=length(thetahat)
  z=matrix(nrow=(n-r), ncol=m, data=a+rexp(m*(n-r),1/thetahat[iter]))
  simcomplete=rbind(z, matrix(nrow=r, ncol=m, data=rep(yc,m)))
  thetahat=c(thetahat,mean(apply(simcomplete,2,mean)))
  diff=abs(thetahat[iter+1]-thetahat[iter])
}
plot(thetahat, type="l")
thetaMLE=ycbar+(n-r)*a/r
abline(h=thetaMLE, col=2)

### Monte Carlo EM with importance sampling 
thetahat=ycbar
m=100
diff=1
z0=matrix(nrow=(n-r), ncol=m, data=a+rexp(m*(n-r),1/thetahat))
zx0=rbind(z, matrix(nrow=r, ncol=m, data=rep(yc,m)))
w=matrix(nrow=n, ncol=m)
while(diff>10^-6)
{
  iter=length(thetahat)
  for(i in 1:m){w[,i]=prod(dexp(zx0[,i],1/thetahat[iter]))/prod(dexp(zx0[,i],1/thetahat[1]))
                w[,i]=n*w[,i]/sum(w[,i])}
  zj=w*zx0
  thetahat=c(thetahat,mean(apply(zj,2,mean)))
  diff=abs(thetahat[iter+1]-thetahat[iter])
}
plot(thetahat, type="l")
thetaMLE=ycbar+(n-r)*a/r
abline(h=thetaMLE, col=2)