import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def MH(x0, f, dprop, rprop, N, B, x):
    """
    ## x0: initializing values for MH algorithm
    ## f: the pdf of posterior
    ## dprop: the pdf of the proposal distribution
    ## rprop: the random number gen for proposal
    ## N: Number of samples needed
    ## B: Number of burn in samples
    ## function generates N+B samples then drops B
    """
    x = np.zeros((N+B, len(x0)))
    fx = np.zeros(N+B)
    x[0,:] = x0
    fx[0] = f(x0)
    counter = 0

    for i in range(1, (N+B)):
        u = rprop(x[i-1,])
        fu = dpost(u)
        r = np.log(fu) + np.log(dprop(x[i-1,])) - np.log(fx[i-1]) - np.log(dprop(u))
        #r <- log(fu) + log(dprop(x[i-1,])) - log(fx[i-1]) - log(dprop(u))

        Rho = min([np.exp(r), 1])
        unif = np.random.uniform(0,1,1)
        if unif <= Rho:
            counter += 1 ## update counter
            x[i,] = u
            fx[i] = fu
        else:
            x[i,] = x[i-1,]
            fx[i] = fx[i-1]
    results = {
        'x':x[(B+1):,], 
        'fx':fx[(B+1):], 
        'rate': counter / (N + B)
        }
    return results



# Posterior distribution
def dpost(theta):
    a = theta[0]
    b = theta[1]
    p = 1 - 1 / (1 + np.exp(a + b * x)) ## logistic CDF
    lik = np.exp(sum(dbinom(y, size=1, prob=p, log=True)))
    dprior = np.exp(a) * np.exp(-np.exp(a) / b_mme) * 1/ b_mme
    return (lik * dprior)



## normal pdf
def dnorm(x, mean, sd):
    out = stats.norm.pdf(x, loc=mean, scale=sd)
    return out

## binomial pdf
def dbinom(x, size, prob, log=True):
    out = stats.binom.logpmf(k=x, n=size, p=prob)
    return out


## density of proposal
def dprop(theta):
    """
    Proposal distribution (independent proposal, so "theta0" is not used)
    """
    a = theta[0]
    b = theta[1]
    pr1 = np.exp(a) * np.exp(-np.exp(a) / b_mme) / b_mme
    pr2 = dnorm(b, b_mle, np.sqrt(var_b_mle))
    return (pr1 * pr2)




def rprop(theta0):
    a = np.log(np.random.exponential(1/b_mme))
    b = np.random.normal(b_mle, np.sqrt(var_b_mle), 1)
    return (a, b)



if __name__ == "__main__":
    # test Data
    y = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0])
    x = np.array([53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70, 70, 70, 72, 73, 75, 75,
        76, 76, 78, 79, 81])

    ## starting values
    a_mle = 15
    b_mle = 0.25
    var_a_mle = 7.37
    var_b_mle = 0.1082
    b_mme = np.exp(a_mle + 0.577216)



    # Run Metropolis-Hastings
    N = 30000 ## Final sample size
    B = 5000 ## burn-in
    mh_out = MH((a_mle, b_mle), dpost, dprop, rprop, N, B, x = x)
    alphaMH = mh_out['x'][:,0]
    betaMH = mh_out['x'][:,1]

    print("""
        MAP parameter estimates:
        alpha: %f
        beta: %f
        """ % (np.mean(alphaMH), np.mean(betaMH)))


    Routput = """
    Coefficients:
                Estimate Std. Error z value Pr(>|z|)  
    (Intercept)  15.0429     7.3786   2.039   0.0415 *
    temperature  -0.2322     0.1082  -2.145   0.0320 *
    """


    print("Compared with output from R: %s" % (Routput))

    plt.hist(alphaMH, 50, normed=1, facecolor='gray', alpha=0.75)
    plt.title("alphaMH Histogram")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.show()


    plt.hist(betaMH, 50, normed=1, facecolor='gray', alpha=0.75)
    plt.title("BetaMH Histogram")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.show()
