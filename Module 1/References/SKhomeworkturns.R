tryrolls=function(N)
{
  count=rep(0,N)
for(i in 1:N){
  x=rep(0,11)
while(sum(x)<11)
{
count[i]=count[i]+1
roll=sum(sample(1:6,2,replace=T))
x[roll-1]=1
}
}
return(count)
}
output=tryrolls(10000)
mean(output)
hist(output, plot=T, breaks=c(min(output):max(output)))