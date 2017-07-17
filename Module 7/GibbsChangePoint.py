import numpy as np
import matplotlib.pyplot as plt


disasters_array = np.array([4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6,
                         3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
                         2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0,
                         1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
                         0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2,
                         3, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
                         0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1])



plt.figure(figsize=(12.5, 3.5))
n_count_data = len(disasters_array)
plt.bar(np.arange(1851, 1962), disasters_array, color="#348ABD")
plt.xlabel("Year")
plt.ylabel("Disasters")
plt.title("UK coal mining disasters, 1851-1962")
plt.xlim(1851, 1962);
plt.show()


## Poisson plots
fig, axes = plt.subplots(1, 4, figsize=(18,4), sharex=True, sharey=True)
for i,l in enumerate([1, 5, 12, 25]):
    axes[i].hist(np.random.poisson(l, 1000), histtype="stepfilled")
    axes[i].annotate(r'$\lambda$=%i' % l, xy=(1, 1), xytext=(30, 350), fontsize=20)




## gamma plots

fig, axes = plt.subplots(1, 4, figsize=(18,4))
for i,p in enumerate([(1, 10), (1, 100), (10, 10), (0.1, 100)]):
    axes[i].hist(np.random.gamma(*p, size=1000), histtype="stepfilled")
    axes[i].set_xlabel(r'$\alpha$=%i, $\beta$=%i' % (p[0], p[1]), fontsize=18)



# Function to draw random gamma variate
rgamma = np.random.gamma

def rcategorical(probs, n=None):
    # Function to draw random categorical variate
    return np.array(probs).cumsum().searchsorted(np.random.sample(n))




def dgamma(lam, a, b):
    return lam**(a-1) * np.exp(-b*lam)


alpha = 1.0
beta = 10


# Specify number of iterations
n_iterations = 10**4

# Initialize trace of samples
mu, lambda2, tau = np.zeros((3, n_iterations+1), dtype=int)


mu[0] = 6
lambda2[0] = 2
tau[0] = 50



# Sample from conditionals
for i in range(n_iterations):
    
    # Sample early mean
    if tau[i] > 110:
        mu[i+1] = rgamma(disasters_array.sum() + alpha, 1./(tau[i] + beta))
    else:
        mu[i+1] = rgamma(disasters_array[:tau[i]].sum() + alpha, 1./(tau[i] + beta))

    # Sample late mean
    lambda2[i+1] = rgamma(disasters_array[tau[i]:].sum() + alpha, 
                          1./(n_count_data - tau[i] + beta))
    if lambda2[i+1] < 1:
        lambda2[i+1] = rgamma(disasters_array.sum() + alpha, 
                          1./(n_count_data - tau[i] + beta))
    
    # Sample changepoint: first calculate probabilities (conditional)
    #p = np.array([dgamma(lambda1[i+1], disasters_array[:t].sum() + alpha, t + beta)*
    #         dgamma(lambda2[i+1], 
    #                disasters_array[t:].sum() + alpha, 
    #                n_count_data - t + beta)
    #         for t in range(n_count_data)])
    # ... then draw sample

    p = np.array(
        [ np.exp( (lambda2[i+1] - mu[i+1]) * t)  * (mu[i+1]/lambda2[i+1])**np.sum(disasters_array[:t])   for t in range(n_count_data) ]
    )
    tau[i+1] = rcategorical(p/p.sum())
    #tau[i+1] = np.random.randint(dtype=int, low=, high=111)

    if i < 5:
        print mu[:i+1]
        print lambda2[:i+1]
        print tau[:i+1]
        print "p"
        print p
        




param_dict = {r'$\lambda_1$':mu, r'$\lambda_2$':lambda2, r'$\tau$':tau}
for p in param_dict:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(param_dict[p][100:])
    axes[0].set_ylabel(p, fontsize=20, rotation=0)
    axes[1].hist(param_dict[p][n_iterations/2:])

plt.show()








