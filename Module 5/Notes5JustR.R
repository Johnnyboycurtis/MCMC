### R code for Module 5 examples

##Numerical optimization using the optimize funtion in R

  xm=rcauchy(500)  ## random draws from Cauchy
f=function(theta){-sum(log(1+(x-theta)^2))}  #f is the log likelihood
## We will maximize f with respect to theta given our observed xs
thetamle=c()
for (i in 1:500)  ## Compute the sequence of maximizations based on increasing sample size
{
  x=xm[1:i]
  thetamle[i]=optimize(f,interval=c(-10,10), maximum=T)$max
}
plot(thetamle, type="l")
abline(h=0, col=2, lty=2)

  
#Here is what happens when we use the likelihood instead of the log likelihood.  The likelihood becomes too small as n increases and the estimate of the MLE diverges.
  
  
    xm=rcauchy(500)  ## random draws from Cauchy
  f=function(theta){prod(dcauchy(x,theta,1))}  #f is the likelihood
  ## We will maximize f with respect to theta given our observed xs
  thetamle=c()
  for (i in 1:500)  ## Compute the sequence of maximizations based on increasing sample size
  {
    x=xm[1:i]
    thetamle[i]=optimize(f,interval=c(-10,10), maximum=T)$max
  }
  plot(thetamle, type="l")
  abline(h=0, col=2, lty=2)
  
  
  ### Newton-Raphson optimization using nlm in R for mean and variance of normal
  
  xvec = c(2,5,3,7,-3,-2,0) # or some other numbers
  #then define a function (which is negative of the log lik)
  
  fn = function(theta) {
    sum ( 0.5*(xvec - theta[1])^2/theta[2] + 0.5* log(theta[2]) )
  }
  nlm(fn, p= c(0,1), print.level=1)
  mean(xvec)
  var(xvec)*(length(xvec)-1)/length(xvec)
  mmt=c()
  plot(mean(xvec), var(xvec),xlim=c(-2,12), ylim=c(5,50))
  for(i in 1:(nlm(fn,c(0,1))$it))
  {
    mmt=rbind(mmt, nlm(fn, c(0,1), iterlim=i)$est)
  }
  lines(mmt, lwd=2)
  mmt
  
  ## Stochastic search to find max of a trig function on interval 0,1
  f=function(x){(cos(50*x)+sin(20*x))^2}
  curve(f,0,1)
  Nsim=10^3
  u=runif(Nsim)
  y=f(u)
  xstar=which.max(y)
  print(c(u[xstar], f(u[xstar])))
  fstar=cummax(y)
  plot(fstar, type="l", ylim=c(3,3.85))
  cumWhich <- function (x, FUN = which.max) { 
    sapply(seq_len(length(x)), function (i) FUN(x[1:i]))}
  ustar=u[cumWhich(y)]
  plot(ustar,pch=16)  
  
  
  ##illustration of temperature in Boltzman-Gibbs transforms for mixture of betas
  f=function(x){.3*dbeta(x,20,80)+.5*dbeta(x,50,50)+.2*dbeta(x,80,20)}  #mixture of betas
  expf=function(x, te){exp((.3*dbeta(x,20,80)+.5*dbeta(x,50,50)+.2*dbeta(x,80,20))/te)} # density in exponential
  expft=function(x, te){expf(x, te)/integrate(expf,0,1,te=te)$value} #normalized value
  myseq=seq(0,1,length=501)
  curve(f,  ylim=c(0,15))
  lines(myseq, expft(myseq, te=2), col=2)
  lines(myseq, expft(myseq, te=3), col=3)
  lines(myseq, expft(myseq, te=.8), col=4)
  lines(myseq, expft(myseq, te=.6), col=5)
  legend(x="topleft", lty=1, lwd=2, col=c(1:6), legend=c("1", "2", "3", ".8", ".6"), bty="n")
  
  ## Simulated annealing example with mixture of betas
  h=function(x){.3*dbeta(x,20,80)+.5*dbeta(x,50,50)+.2*dbeta(x,80,20)}  #mixture of betas
  
  cool1=function(iter){1/log(1+iter)}
  cool2=function(iter){1/(1+iter)^2}
  
  scale1=function(iter){sqrt(cool1(iter))}
  scale2=function(iter){sqrt(cool2(iter))}
  
  ##Do the annealing as a function.  x is the starting value, theta_0. cool is 
  ##a function with the cooling schedule and scale is a function with the scale parameter for the uniform proposal density
  SAbeta=function(x=runif(1),cool=cool1, sc=scale1)
  {
    
    hval=hcur=h(x) 
    diff=iter=1
    temp=cool(iter)
    scale=sc(iter)
    
    while(diff>10^(-4)){
      prop=x[iter]+runif(1,-1,1)*scale  ##proposed new value of theta
      if((prop>1)||(prop<0)||(log(runif(1))*temp>(h(prop)-hcur)))prop=x[iter] ##reasons to repeat theta_(t-1) instead of taking new prop
      x=c(x,prop)
      hcur=h(prop)
      hval=c(hval, hcur)
      if((iter>100))diff=max(hval)-max(hval[1:(iter/2)])
      temp=cool(iter)
      scale=sc(iter)
      iter=iter+1
    }
    return(cbind(x,hval))
  }
  
  
  
  curve(h, lwd=2)
  output=SAbeta(x=runif(1),cool=cool1, sc=scale1)
  lines(output[,1],output[,2])
  maxspot=which.max(output[,2])
  #print(output)
  print(output[maxspot,])
  
  
  ###Experiment with different scales/cooling to see if they end up in the right spot.
  
  rightmax=function(cool,scale)
  {
    output=SAbeta(x=runif(1),cool, sc=scale)
    return(round(output[which.max(output[,2]),1],1)==.5)
  }
  c1s1results=replicate(100,rightmax(cool1,scale1))
  mean(c1s1results)
  c2s2results=replicate(100,rightmax(cool2,scale2))
  mean(c2s2results)
  scale3=function(iter){3*sqrt(cool2(iter))}
  c2s3results=replicate(100,rightmax(cool2,scale3))
  mean(c2s3results)
  
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
  r=sum(x<a)
  yc=y[1:r]
  ycbar=mean(yc)
  ### EM
  thetahat=cur=ycbar
  diff=1
  while(diff>10^-4)
  {
    thetahat=c(thetahat,(ycbar*r+(n-r)*(a+thetahat[length(thetahat)]))/n )
    diff=abs(thetahat[length(thetahat)]-thetahat[length(thetahat)-1])
  }
  plot(thetahat, type="l")
  thetaMLE=ycbar+(n-r)*a/r
  abline(h=thetaMLE, col=2)
  
  
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