MCPiEstimate=function(N=1000, pl=TRUE) { ## N is the number of random draws to take inside the unit square.  If no N is specified, the default is 1000.  pl is an indicator of whether to plot the points
x=runif(N) ## draw the x coordinate randomly and uniformly from the interval (0,1)
y=runif(N) ## draw the y coordinate randomly and uniformly from the interval (0,1)
incircle=(sqrt(x^2+y^2)<1) ## Determine whether each point is in the unit circle
par(pty="s") ## To tell the plot region to be square - it is rectangular by default
if(pl){plot(x,y,col=incircle+1,pch=16)} #plot the x,y pairs with incirle in red (col=2) and out of cirle in black (col=1)
p=sum(incircle/N) ## proportion of random points in the circle
return(p*4)} ## Monte Carlo estimate of pi, output of the function

MCPiEstimate(1000)
MCPiEstimate(10000)
MCPiEstimate(100000)
pi1000=replicate(500,MCPiEstimate(1000, pl=F))
pi10000=replicate(500,MCPiEstimate(10000, pl=F))
pi1000
pi10000
mean(pi1000)
mean(pi10000)
sd(pi1000)
sd(pi10000)
plot(density(pi10000),xlim=c(3,3.3))
lines(density(pi1000))