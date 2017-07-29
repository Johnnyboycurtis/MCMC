


### Bayesian Analysis with Normal Data

Suppose we have a random sample $y_1, ..., y_n | \mu, \tau \sim N(\mu, 1/\tau)$. Here we use a precision parameter $\tau$ instead of $\sigma$ for standard devation (or $\sigma^2$ for variance) .
Place two independent priors on the unkown parameters $\mu$ and $\tau$.

$$
\mu \sim N(a, 1/b) \perp \tau \sim Gamma(c,d)
$$

The posterior is given by
$$
p(\mu, \tau | \mathbf{y} ) = \frac{L(\mu, \tau | \mathbf{y}) \times p(\mu, \tau) }{m(\mathbf{y})}
$$

where $p(\mu, \tau) = p(\mu) \times p(\tau)$ and $m(\mathbf{y}) = \int L(\mu, \tau | \mathbf{y}) \times p(\mu, \tau)$
  

The likelihood is 
$$
  \begin{aligned}
L(\mu, \tau | \mathbf{y} ) & = \prod_{i=1}^{n} f(y_i | \mu, \tau) \\
\ & = \frac{\tau}{2\pi}^{n/2} exp(\frac{\tau}{2} \sum (y_i - \mu)^2 ) \\
\ & = \frac{\tau}{2\pi}^{n/2} exp{ \frac{\tau}{2} (n-1)s^2 - \frac{\tau}{2}(\bar{y} - \mu)^2 } \text{ ( obtained after some algebra )}  \\ 
\ & \propto \tau^{n/2} exp{ \frac{\tau}{2} (n-1)s^2 - \frac{\tau}{2}(\bar{y} - \mu)^2 } \text{ ( drop unimportant constants ) }  \\
\end{aligned}
$$
  


The posterior is then given by
$$
\begin{aligned}
p(\mu, \tau) & \propto  L(\mu, \tau | \mathbf{y}) \times p(\mu, \tau) \\
\ & = \tau^{n/2} exp( \frac{\tau}{2} (n-1)s^2 - \frac{\tau}{2} (\bar{y} - \mu)^2 ) \times \frac{d^c}{ \Gamma(c) } \tau^{c-1} exp(-b \tau) \times \sqrt{\frac{b}{2 \pi}} exp(- \frac{b}{2} (\mu - a)^2 ) \\
\ & \text{ ( next, drop unimportant constants ) } \\
\ & \propto \tau^{n/2} exp( \frac{\tau}{2} (n-1)s^2 - \frac{\tau}{2} (\bar{y} - \mu)^2 ) \times \tau^{c-1} exp(-b \tau) \times \ exp(- \frac{b}{2} (\mu - a)^2 ) \\
\ & = \tau^{n/2 + c - 1} exp( \frac{\tau}{2} (n-1)s^2 - \frac{\tau}{2} (\bar{y} - \mu)^2 ) \times exp(-b \tau) \times \ exp(- \frac{b}{2} (\mu - a)^2 ) \\
\end{aligned}
$$
  
And the final result is some function that is not recognizable as any commonly known distribution. Hence, we're stuck if we are to try to sample from this distribution directly!

However, if we are too look at the conditional posteriors, $\mu | \tau, \mathbf{y}$ 
$$
p(\mu | \tau, \mathbf{y}) \propto exp(- \frac{1}{2} (n \tau + b) (\mu - \hat{\mu}_{\tau} )^2 )
$$
and $\tau | \mu,\mathbf{y}$ 
$$
p(\tau | \mu, \mathbf{y}) \propto \tau^{(n/2 + c) - 1} exp(- \tau (d + \frac{n (\bar{y} - \mu)^2 + (n-1)s^2 }{2} ))
$$
we find that there are recognizable distributions. They are $N(\hat\mu, \frac{1}{n \tau + b})$ and $Gamma(\frac{n}{2} + c, d +  \frac{n (\bar{y} - \mu)^2 + (n-1)s^2 }{2})$

Given that we can sample from the conditional distributions, we can utilize a Gibbs sampler to sample from the posterior distribution $p(\mu, \tau | \mathbf{y})$.






```{r}

y = mtcars$mpg

plot(density(y), main = "Density of MPG", 
     lty = 2, xlab = "MPG")



# This function takes the prior mean and a percentile and returns the corresponding hyperparameter values for the normal prior
find.normal <- function(prior.mean, percentile, p=0.95){
  prior.sd <- (qnorm(p) / (percentile - prior.mean))^{-1}
  return(list(mean=prior.mean, sd=prior.sd, precision=1/prior.sd^2)) 
}


# This function takes the prior mode and a percentile of sigma and returns the corresponding parameters of the Gamma distribution for tau 
find.tau.gamma <- function(prior.sigma.mode, sigma.percentile, p=0.95, start.a=1, end.a=1000, increment=0.01){
  a <- seq(from=start.a, to=end.a,  by=increment)
  b <- prior.sigma.mode^2 * (a + 1)
  i <- which.max(qgamma(1-p, a, b) > 1/sigma.percentile^2)
  return(list(a=a[i], b=b[i])) 
}


```


```{r}

N = 3 * 10^3
n = length(y)
ybar = mean(y)
s = sd(y)
mu = numeric(N)
tau = numeric(N)

## initial values
mu = ybar
tau = 1/sd(y)


## Normal hyper parameters
a <- find.normal(prior.mean=21, percentile = 24, p=0.95)$mean
b <- find.normal(prior.mean=21, percentile = 24, p=0.95)$precision

## Gamma Hyper parameters
gamma_HyperParams <- find.tau.gamma(prior.sigma.mode= 1.2, sigma.percentile = 2.462, p=0.95) # Returns shape and rate parameters for the Gamma distribution of tau
gamma_HyperParams

c = gamma_HyperParams$a ## arbitrary
d = gamma_HyperParams$b ## arbitrary

for(i in 2:N){
  SHAPE = c + n/2
  RATE = d + 0.5*( n*(ybar - mu[i-1])^2 + (n-1)*s^2  )
  tau[i]  = rgamma(n = 1, shape = SHAPE, rate = RATE)
  
  MU = mu[i-1]
  SD = 1/(n*tau[i] + b)
  mu[i] = rnorm(n = 1, mean = MU, sd = SD)  
}

B = 10^2

muPosterior = mu[-(1:B)]
tauPosterior = tau[-(1:B)]

plot(density(mu), main = "Density of mu", 
     lty = 2, xlab = "mu")




plot(density(tau), main = "Density of tau", 
     lty = 2, xlab = "tau")



```






