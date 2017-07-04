# Ex 1.8
set.seed(829125)
x=seq(-3,3,le=5)
y=2+4*x+rnorm(5)
lm(y~x)
fit=lm(y~x)
Resdata=fit$residuals
nBoot=2000
B=array(0,dim=c(nBoot,2))
for(i in 1:nBoot){
  ystar=y+sample(Resdata, replace=T)
  Bfit=lm(ystar~x)
  B[i,]=Bfit$coefficients
}
Bciint<-quantile(B[,1],c(.025, .975))  #Bootstrap confidence interval for intercept
Bcislope<-quantile(B[,2], c(.025, .975))
rbind(Bciint, Bcislope)
confint(fit) #t distribution confidence intervals