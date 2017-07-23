
study = matrix(data = c(21, 2,
                        15, 3), nrow = 2, ncol = 2, byrow = TRUE,
               dimnames = list(c("surgery", "radiation"), 
                               c("controlled", "not controlled")))
alpha = 0.05
B = 10^4


## SOME SET UP
x = study
almost.1 <- 1 - 64 * .Machine$double.eps


n <- sum(study)
sr <- rowSums(x)
sc <- colSums(x)
E <- outer(sr, sc, "*")/n ## ORDER MATTERS
dimnames(E) <- dimnames(study)


simulateChisq <- function(B, E, sr, sc){
    results = numeric(B)
    for(i in 1:B){
        dat = unlist(r2dtable(1, sr, sc))
        M = matrix(dat, ncol = length(sc), nrow = length(sr))
        val = sum( sort( (M - E)^2 / E, decreasing = TRUE))
        results[i] = val
    }
    return(results)
}

#setMETH()
#tmp <- .Call(C_chisq_sim, sr, sc, B, E)
#tmp <- .Call(stats:::C_chisq_sim, sr, sc, B, E, PACKAGE="stats")
tmp <- simulateChisq(B, E, sr, sc)
STATISTIC <- sum(sort((x - E)^2/E, decreasing = TRUE))
PARAMETER <- NA
PVAL <- (1 + sum(tmp >= almost.1 * STATISTIC))/(B + 1)


ChisqTest <- function(data, Simulations){
    ## data should be a 2X2 matrix
    x = data
    B = Simulations
    #almost.1 <- 1 - 64 * .Machine$double.eps
    n <- sum(x)
    sr <- rowSums(x)
    sc <- colSums(x)
    E <- outer(sr, sc, "*")/n ## ORDER MATTERS
    dimnames(E) <- dimnames(study)
    tmp <- simulateChisq(B, E, sr, sc)
    Stat <- sum(sort((x - E)^2/E, decreasing = TRUE))
    pval <- (1 + sum(tmp >= almost.1 * STATISTIC))/(B + 1)
    df = 2 ## only option for this example
    rawPVal = pchisq(q = Stat, df = df, lower.tail = FALSE)
    out = list(PearsonStat = Stat, MonteCarloPVal = pval, rawPVal = rawPVal)
    return(out)
}


plot(density(tmp))
print(STATISTIC)
print(PVAL)


w = chisq.test(study, simulate.p.value=TRUE, B = B)
print(w)






library(DescTools)

G_test_stat = GTest(study)
print("GTest:")
print(G_test_stat)

GStat <- function(x){
  ## x should be a matrix
  n = sum(x)
  O = x
  SC = colSums(O)
  SR = rowSums(O)
  E = outer(SR, SC, "*") / n ## order of sr, sc matters
  G = 2*( (O[1,1]*log(O[1,1]/E[1,1])) + 
            (O[1,2]*log(O[1,2]/E[1,2])) + 
            (O[2,1]*log(O[2,1]/E[2,1])) + 
            (O[2,2]*log(O[2,2]/E[2,2]))  )
  #G = GTest(m)$statistic
  return(G)
}

print("GStat: ")
print(GStat(study))



simulateLogLike <- function(B, E, sr, sc){
    results = numeric(B)
    for(i in 1:B){
        dat = unlist(r2dtable(1, sr, sc))
        O = matrix(dat, ncol = length(sc), nrow = length(sr))
        G = 2*( (O[1,1]*log(O[1,1]/E[1,1])) + 
            (O[1,2]*log(O[1,2]/E[1,2])) + 
            (O[2,1]*log(O[2,1]/E[2,1])) + 
            (O[2,2]*log(O[2,2]/E[2,2]))  )
        results[i] = G
    }
    return(results)
}

