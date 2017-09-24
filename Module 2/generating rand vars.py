import numpy as np
import matplotlib.pyplot as plt

## set style
plt.style.use("ggplot") ## use plt.style.availalbe to see all styles


u = np.random.uniform(low=0,high=1,size=2000)
plt.hist(u, bins = 30, normed=1, color='blue', alpha = 0.5)
plt.title("Uniform samples")
plt.xlabel("Uniform samples")
plt.ylabel("Probability")
plt.show()


import scipy.stats as stats

stats.probplot(u, dist="uniform", plot=plt)
plt.ylabel("Ordered Sample Values")
plt.title("Q-Q Plot")
plt.show()



# import statsmodels.api as sm
# test = np.random.normal(0,1, 1000)
# sm.qqplot(test, line='45')
# plt.show()


#acorr(x, normed=True, detrend=mlab.detrend_none, usevlines=True,
#      maxlags=10, **kwargs)

plt.acorr(u, normed=True, maxlags = 10)
plt.title("ACF Plot")
plt.xlabel("Lags")
plt.ylabel("autocorrelation")
plt.show()



from statsmodels.graphics.tsaplots import plot_acf
#1. pandas
plot_acf(u)
plt.title("ACF Plot")
plt.xlabel("Lags")
plt.ylabel("autocorrelation")
plt.show()



# u = runif(2000) ## generate a sample of unif(0,1) samples

# hist(u, probability=TRUE, col = "gray", border = "white") ## plot histogram
# ## Q-Q plot for `runif` data against true theoretical distribution:
# qqplot(x = qunif(ppoints(500)), y = u,  main = expression("Q-Q plot for Unif(0,1)"))
# qqline(y = u, distribution = qunif, prob = c(0.1, 0.6), col = 2)

# acf(x = u, main = "autocorrelation of samples")  ## autocorrelation function


stats.probplot(u, dist="uniform", plot=plt)
plt.title("(1-U) Uniform samples")
plt.xlabel("(1-U) Uniform samples")
plt.ylabel("Probability")
plt.show()



# qqplot(qunif(ppoints(500)),1-u, main = "U against 1 - U")
# qqline(u, distribution = qunif, col = 2)


v = 10 * np.random.uniform(low=0.0, high=1.0, size = 1000)
plt.hist(v, bins=20, color = "gray", alpha = 0.5, normed=True)
plt.title("V = 10 x U samples")
plt.xlabel("V")
plt.ylabel("Density")


w = np.random.uniform(low=0.0, high=10.0, size = 1000)
plt.hist(w, bins=20, color = "green", alpha = 0.5, normed=True)
plt.title("W = Unif(0, 10) samples")
plt.xlabel("W")
plt.ylabel("Density")


# v = 10 * runif(1000, 0, 1) ## v ~ unif(0, 10)
# hist(v, col = "gray", border = "white")
# w = runif(1000, 0, 10) ## unif(0, 10)
# hist(w, col = "gray", border = "white")



