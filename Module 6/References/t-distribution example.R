#!/usr/bin/R

## Student's t density with v degres of freedom is given by *pdf*

#Calculate the mean of a t distribution with v = 4 degrees of freedom using a M-H algorithm with 
#candidate densities N(0,1) and t_{v = 2}

set.seed(123)
N = 10^6

#dt(x = x, df = 4)

X = numeric(N)
X[1] = rnorm(1) ## initialize the starting value

for(i in 1:N){
    Y = rnorm(1) ## independent of X_i
    rho = (dt(Y, df = 4) * dnorm(X[i])) /
            (dt(X[i], df = 4) * dnorm(Y))
    U = runif(1)
    if(U <= rho){
        X[i+1] = Y
    } else{
        X[i+1] = X[i]
    }

}




plot(density(X), type = "l", 
        lty = 2, main = "M-H with N(0,1) candidate")
curve(dt(x, df = 4), add = TRUE, col = 4)



## now with a t distribution with df = 2

X = numeric(N)
X[1] = rt(1, df = 2) ## initialize the starting value

for(i in 1:N){
    Y = rt(1, df = 2) ## independent of X_i
    rho = (dt(Y, df = 4) * dt(X[i], df = 2)) /
            (dt(X[i], df = 4) * dt(Y, df = 2))
    U = runif(1)
    if(U <= rho){
        X[i+1] = Y
    } else{
        X[i+1] = X[i]
    }

}




plot(density(X), type = "l", 
        lty = 2, main = "M-H with t_df = 2 candidate")
curve(dt(x, df = 4), add = TRUE, col = 4)
