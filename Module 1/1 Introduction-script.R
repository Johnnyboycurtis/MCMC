
study = matrix(data = c(21, 2,
                        15, 3), nrow = 2, ncol = 2, byrow = TRUE,
               dimnames = list(c("surgery", "radiation"), 
                               c("controlled", "not controlled")))

print(study)
## set up

## function will generate chi-squared statistics 
## using the expected distribution of the data
simulateChisq <- function(B, E, sr, sc){
    results = numeric(B)
    for(i in 1:B){
      ## review r2dtable documentation
        dat = unlist(r2dtable(1, sr, sc)) ## simulated contingency table
        M = matrix(dat, ncol = length(sc), nrow = length(sr))
        val = sum( sort( (M - E)^2 / E, decreasing = TRUE))
        results[i] = val
    }
    return(results)
}


#[`r2dtable`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/r2dtable.html)

#[Random generation of a table ](http://people.sc.fsu.edu/~jburkardt/f_src/asa159/asa159.html)


## Monte Carlo Hypothesis Testing, pt. 5 


ChisqTest <- function(data, Simulations){
    x = data ## data
    B = Simulations ## number of simulations to generate
    n <- sum(x) ## total number of observations
    sr <- rowSums(x) ## sum of rows
    sc <- colSums(x) ## sum of cols
    E <- outer(sr, sc, "*")/n ## ORDER MATTERS
    dimnames(E) <- dimnames(study)
    tmp <- simulateChisq(B, E, sr, sc) ## simulated data
    Stat <- sum(sort((x - E)^2/E, decreasing = TRUE)) ## chi^2 statistic
    pval <- (1 + sum(tmp >=  Stat))/(B + 1) ## MC pvalue
    rawPVal = pchisq(q = Stat, df = 2, lower.tail = FALSE)
    out = list(PearsonStat = Stat, 
               MonteCarloPVal = pval, 
               rawPVal = rawPVal)
    return(out)
}





#We then generate our test statistics.


#set.seed(123)

results <- ChisqTest(study, 10000)

print(results)

## compare against chisq.test()






## Bayesian Example

N = 10^4
set.seed(123)
x = rbeta(n = N, shape1 = 576 + 1, shape2 = 1028 - 576 + 1)
d = density(x)
hist(x = x, probability = TRUE, 
      main = "Beta Posterior Distribution",
     xlab = expression(theta), ylab = "Density",
     ylim = c(0,40), col = "gray", border = "white")
lines(x = d$x , y = d$y, type = "l", col = 2)
abline(v = median(x), lty = 3, col = "3")

print("Median: ")
print(quantile(x = x, probs = c(0.025, 0.5, 0.975)))




