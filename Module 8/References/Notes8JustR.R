library(mcsm)
randogibsk=function (T = 10^3) 
{
  ###Simulate data from the random effects logit model.
  beta0 = -3
  sigma0 = 1
  n = 20   ## number of groups
  m = 35   ## number of observations per group
  x = matrix(sample(c(-1, 0, 1), n * m, rep = T), nrow = n)  ## covariate value
  y = matrix(0, ncol = m, nrow = n)
  tru = rnorm(n)   ### true random effect values
  for (i in 1:n) y[i, ] = (runif(m) < 1/(1 + exp(-beta0 * x[i,   ## observed y values
                                                            ] - tru[i])))
  mlan = as.numeric(glm(as.vector(y) ~ as.vector(x) - 1, family = binomial)$coef)  ##maximum likelihood est of beta
  
  ## functions for the MH algorithm
  likecomp = function(beta, sigma, u) {
    sum(as.vector(y) * (beta * as.vector(x) + rep(u, m) - 
                          rep(tru, m))) - sum(log(1 + exp(beta * as.vector(x) + 
                                                            rep(u, m)))) - sum(u^2)/(2 * sigma^2) - log(sigma)
  }
  gu = function(mu, i, beta, sigma) {           ### likelihood for a single group
    sum((y[i, ] * (beta * x[i, ] + mu)) - log(1 + exp(beta * 
                                                        x[i, ] + mu))) - 0.5 * mu^2/sigma^2
  }
  pro = function(beta, u) {
    exp(as.vector(y) * (beta * as.vector(x) + rep(u, m)))/(1 + 
                                                             exp(as.vector(y) * (beta * as.vector(x) + rep(u, 
                                                                                                           m))))^2
  }
  
  ##initialize and storage
  beta = mlan
  sigma = factor = 1
  acpt = bcpt = 0
  u = rnorm(n)
  samplu = matrix(u, nrow = n)
  for (iter in 2:T) {
    ##MH for random effects
    u = rnorm(n)
    for (i in 1:n) {
      mu = samplu[i, iter - 1]
      u[i] = factor * sigma[iter - 1] * rnorm(1) + mu
      if (log(runif(1)) > gu(u[i], i, beta[iter - 1], sigma[iter - 
                                                              1]) - gu(mu, i, beta[iter - 1], sigma[iter - 
                                                                                                      1])) {
        u[i] = mu
      }
      else {
        acpt = acpt + 1
      }
    }
    samplu = cbind(samplu, u)
    ## update variance
    sigma = c(sigma, 1/sqrt(2 * rgamma(1, 0.5 * n)/sum(u^2)))
    tau = sigma[iter]/sqrt(sum(as.vector(x^2) * pro(beta[iter - 
                                                           1], u)))
    ## MH for coefficient
    betaprop = beta[iter - 1] + rnorm(1) * factor * tau
    if (log(runif(1)) > likecomp(betaprop, sigma[iter], u) - 
          likecomp(beta[iter - 1], sigma[iter], u)) {
      betaprop = beta[iter - 1]
    }
    else {
      bcpt = bcpt + 1
    }
    beta = c(beta, betaprop)
    ## adjust scale of proposal distribution
    if (iter > 100) {
      if (bcpt < 0.1 * iter) 
        factor = factor/3
      if (bcpt > 0.9 * iter) 
        factor = factor * 3
    }
  }
  
  ## Check on chain using coda
  library(coda)
  mymcmc=mcmc(cbind(beta, sigma))
  plot(mymcmc)
  cumuplot(mymcmc)
  autocorr.plot(mymcmc)
  list(beta = beta, sigma = sigma)
}
set.seed(1124)
resagibs=randogibsk()

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
    phi<-rgamma(1,c+sum(y[theta.matrix[(i-1),3]:n]),(d+n-theta.matrix[(i-1),3]))
    k<-sample(kseq,1, prob=lapply(kseq,FUN=kprob,y,lambda, phi, a, b, c, d, n))
    theta.matrix<-rbind(theta.matrix, c(lambda, phi, k))
  }
  return(theta.matrix )
}

Nsim=10^4
start<-matrix(c(1,1,40),1,3)
coal.sample<-gibbs.poisson.gamma(start,coal.mining.disasters, Nsim)
#coal.sample<-gibbs.poisson.gamma(start, coal.mining.disasters, 1000)

ksseq=function(beta,M)
{ ks = NULL
  T=length(beta)
  for (t in seq(T/10, T, le = 100)) {
    beta1 = beta[1:(t/2)]
    beta2 = beta[(t/2) + (1:(t/2))]
    beta1 = beta1[seq(1, t/2, by = M)]
    beta2 = beta2[seq(1, t/2, by = M)]
    ks = c(ks, ks.test(beta1, beta2)$p)
  }
  
  plot(seq(1, T, le = 100), ks, pch = 19, cex = 0.7, xlab = "Iterations", 
       ylab = "p-value")
}
par(mfrow=c(1,1))
ksseq(coal.sample[,1],10)

mymcmc=mcmc(coal.sample[1:1000,])
plot(mymcmc)
cumuplot(mymcmc)

Nsim=1500
start2<-matrix(c(.5,.5,20),1,3)
start3<-matrix(c(2,1,55),1,3)
coal.sample1<-gibbs.poisson.gamma(start,coal.mining.disasters, Nsim)
coal.sample2<-gibbs.poisson.gamma(start2,coal.mining.disasters, Nsim)
coal.sample3<-gibbs.poisson.gamma(start3,coal.mining.disasters, Nsim)
coal.sampleparallel=mcmc.list(mcmc(coal.sample1, start=100), mcmc(coal.sample2, start=100), mcmc(coal.sample3, start=100))
gelman.plot(coal.sampleparallel)
gelman.diag(coal.sampleparallel)
acf(coal.sample2)
acfplot(mcmc(coal.sample2))
effectiveSize(mcmc(coal.sample2))

library(coda)
setwd("D://My Documents//STAT 701//Module8Conv")
Dyesvague=read.openbugs("Dyesvague")
Dyesgam=read.openbugs("Dyesgam")
summary(Dyesvague)
plot(Dyesvague, ask=FALSE)
plot(Dyesvague[[1]][1:1000,7],type="l" ,ask=FALSE)
plot(density(Dyesvague[[1]][,3,drop=FALSE]))
cumuplot(Dyesvague[[1]])
autocorr.plot(Dyesvague[[1]],ask=FALSE)
effectiveSize(Dyesvague[[1]])
summary(Dyesgam)
plot(Dyesgam, ask=FALSE)
cumuplot(Dyesgam[[1]])
autocorr.plot(Dyesgam[[1]],ask=FALSE)
effectiveSize(Dyesgam[[1]])
gelman.plot(Dyesvague)
gelman.diag(Dyesgam)
gelman.plot(Dyesgam)
gelman.diag(Dyesgam)