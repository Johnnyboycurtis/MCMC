from scipy.stats import beta
import numpy as np
from numpy import random
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab

## set style
plt.style.use("ggplot") ## use plt.style.availalbe to see all styles

n = 10**4
a=576 + 1
b = 1028-576+1
random.seed(123)
x = random.beta(a=a, b = b, size=n)

# the histogram of the data
n, bins, patches = plt.hist(x, 50, normed=1, facecolor='black', alpha=0.5)

# add a 'best fit' line
y = beta.pdf(bins, a, b)
l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('Posterior samples')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ }\ \beta (a,b)$')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)
median = np.median(x)
plt.vlines(x = median, ymin = 0, ymax = 30, colors='blue', linestyle = 'dashdot')

CI = np.percentile(a = x, q = [0.025, 0.975]) ## 95% CI

print "median: " + str(median)
print "CI: " + str(CI)


plt.show()
