#> source("c:/WorkTransfer/CDCcourse/Labs/MHLogistic.txt",print.eval=TRUE)#
#This analyzes the health policy data#
logdata<-read.table("c:/WorkTransfer/CDCcourse/Labs/LogisticData.txt",header=T)
np<-logdata[,2];met<-logdata[,3];erodd<-logdata[,5]
#Get MLEs, set up densities#
#summary(glm(erodd~met+np, family=binomial(link=logit)))
loglike <- function(a,b,c)sum(erodd*(a+b*met+c*np)-log(1+exp(a+b*met+c*np)))
am<--1.97392;as<- 0.22109;
bm<- 0.28438;bs<-0.09270;
cm<- 0.16221;cs<-0.07971 ;
dcand<-function(a,b,c)dnorm(a,am,as,log=TRUE)+dnorm(b,bm,bs,log=TRUE)+dnorm(c,cm,cs,log=TRUE)
nsim<-1000;ahat<-array(0,dim=c(nsim,1));bhat<-array(0,dim=c(nsim,1));
chat<-array(0,dim=c(nsim,1));
ahat[1]<-rnorm(1,mean=am,sd=as);bhat[1]<-rnorm(1,bm,bs);chat[1]<-rnorm(1,cm,cs)
for (j in 2:nsim)  {
a0<-ahat[j-1];b0<-bhat[j-1];c0<-chat[j-1];
acand<-rnorm(1,mean=am,sd=as);bcand<-rnorm(1,bm,bs);ccand<-rnorm(1,cm,cs)
test<- min(exp(loglike(acand,bcand,ccand)-loglike(a0,b0,c0)+dcand(a0,b0,c0)-dcand(acand,bcand,ccand)),1);
rho<-(runif(1)<test);
ahat[j]<-acand*rho+ahat[j-1]*(1-rho);
bhat[j]<-bcand*rho+bhat[j-1]*(1-rho);
chat[j]<-ccand*rho+chat[j-1]*(1-rho);
}
den<-1:(nsim)
par(mfrow=c(2,3))
hist(ahat,main="Intercept",freq=F,col="green",breaks=30)
par(new=T)
plot(density(ahat),main="",xlab="",ylab="",xaxt="n",yaxt="n")
hist(bhat,main="sickness",freq=F,col="green",breaks=30)
par(new=T)
plot(density(bhat),main="",xlab="",ylab="",xaxt="n",yaxt="n")
hist(chat,main="HMO type",freq=F,col="green",breaks=30)
par(new=T)
plot(density(chat),main="",xlab="",ylab="",xaxt="n",yaxt="n")
amean<-cumsum(ahat)/den;plot(amean,type="l",col="red")
bmean<-cumsum(bhat)/den;plot(bmean,type="l",col="red")
cmean<-cumsum(chat)/den;plot(cmean,type="l",col="red")
