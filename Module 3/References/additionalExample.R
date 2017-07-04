f = function(x){ 
    1.5 * (1 - x^2)
}


curve(f(x), from = 0, to = 1, ylim(0, 2))
myseq = seq(from=0, to=1, by 0.01)
lines(myseq, dbeta(myseq, 1, 2))

Nsim = 10^4
x = rbeta(Nsim, 1, 2)

w = f(x)/dbeta(x, 1, 2)
lines(density(x, weights = w/sum(w), col = 2))
mean(w)
