## script to set environment up


import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


message = """now importing: 
1. numpy as np 
2. matplotlib.pyplot as plt
3. from scipy import stats
"""

print(message)

## set style
plt.style.use("ggplot") ## use plt.style.availalbe to see all styles
plt.style.use('seaborn-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6) ## set plot sizes


def summary(data, quantiles = [0, 25, 50, 75, 100], axis = 0):
    """ To print out summary statistics"""
    titles = " Min.:, 1st Qu.:, Median:, 3rd Qu.:, Max.:".split(",")
    percentiles = np.round(
        np.percentile(a=data, interpolation='midpoint', 
                      q=quantiles, 
                      axis = axis), 4)
    out = dict(zip(titles, percentiles))
    mean = np.round(np.mean(a = data, axis = axis), 4)
    out.update({" Mean:" : mean})
    for k, val in out.items():
        print(k, val)
    #return out

def table(data, prob=False):
    """Count unique values
    return a tuple(value, count)"""
    vals_counts = np.unique(ar=data, return_counts=True)
    if prob:
        density = vals_counts[1] / len(data)
        vals_counts = (vals_counts[0], density)
    out = tuple(zip(*vals_counts))
    return out


def _getcolor(val):
    colors = {0: "#E24A33", 1: "#348ABD", 2: "#988ED5", 3: "#777777", 
              4: "#FBC15E", 5: "#8EBA42", 6: "#FFB5B8", 7: "#E24A33"}
    color_choice = colors.get(val, val)
    return color_choice



def plot(x, y, linestyle='-', color = "#348ABD", title="", xlabel="x", ylabel="y", ylim = None, xlim = None, show=True):
    """Basic line plots"""
    color_choice = _getcolor(color)
    tmp, = plt.plot(x, y, linestyle=linestyle, color=color_choice)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.title(title, fontsize = 14)
    plt.ylabel(ylabel, fontsize = 13)
    plt.xlabel(xlabel, fontsize = 13)
    plt.ylim(ylim)
    plt.xlim(xlim)
    if show:
        plt.show()
        return tmp
    else:
        return tmp
    


def hist(x, bins='doane', xlabel = "X", title = "Histogram",cumulative=False, density=True, color='0.3', edgecolor = 'white', show=True):
    """ Plot Histograms """
    color_choice = _getcolor(color)
    if isinstance(bins, str):
        breaks, _ =  np.histogram(a=x, bins=bins)
        bins = len(bins) * 2

    if density:
        #n = x.shape[0]
        w = np.ones_like(x) #/ n
        normed = True
        ylabel = "Density"
    else:
        w = np.ones_like(x)
        ylabel = "Frequency"

    tmp = plt.hist(x, bins=bins, weights=w, cumulative=cumulative, color=color_choice,  edgecolor=edgecolor, normed=normed)
    plt.xticks(fontsize = 12) 
    plt.yticks(fontsize = 12)
    plt.ylabel(ylabel, fontsize = 13)
    plt.xlabel(xlabel, fontsize = 13)
    plt.title(title, fontsize = 14)
    if show:
        plt.show()
    else:
        return tmp



def g1(x):                                                              
    mu = np.mean(x)                                                     
    sd = np.std(x)                                                      
    return np.mean(((x-mu)/sd)**3)                                      
                                                                         
                                                                         
def sdg1(n):                                                            
    num = 6*(n-2)                                                       
    denom = (n+1)*(n+3)                                                 
    return np.sqrt(num/denom)                                           
                                                                         
                                                                         
def doane(x):                                                          
    n = x.shape[0]                                                     
    g1val = g1(x)                                                      
    sdg1val = sdg1(n)                                                  
    bins = 1 + np.log2(n) + np.log2(1 + np.abs(g1val)/sdg1val)         
    return bins  



def cumMean(data):
    """ calculate cumulative mean """
    m = len(data)
    result = np.cumsum(data) / np.arange(1, m+1)
    return result

def cumSE(data):
    """ calculate cumulative standard error """
    m = len(data)
    mu = np.mean(data)
    num = np.sqrt(np.cumsum((data - mu)**2))
    denom = range(1, m+1)
    result = num/denom ## cummulative mean of (x_i - theta)**2
    return result



