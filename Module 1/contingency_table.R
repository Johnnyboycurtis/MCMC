
study = matrix(data = c(21, 2,
                        15, 3), nrow = 2, ncol = 2, byrow = TRUE,
               dimnames = list(c("surgery", "radiation"), 
                               c("controlled", "not controlled")))
alpha = 0.05
B = 4000

x = study
almost.1 <- 1 - 64 * .Machine$double.eps


n <- sum(study)
sr <- rowSums(x)
sc <- colSums(x)
E <- outer(sr, sc, "*")/n
dimnames(E) <- dimnames(study)


myfun <- function(B, E, sr, sc){
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
tmp <- myfun(B, E, sr, sc)
STATISTIC <- sum(sort((x - E)^2/E, decreasing = TRUE))
PARAMETER <- NA
PVAL <- (1 + sum(tmp >= almost.1 * STATISTIC))/(B + 1)

plot(density(tmp))
print(STATISTIC)
print(PVAL)


w = chisq.test(study, simulate.p.value=TRUE)
print(w)






library(DescTools)

GTest(study)

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


GStat(study)
