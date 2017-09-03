n = 50000 ## 20000 is too small!!
x = runif(n)

h = function(x){
  v = (cos(50*x) + sin(50*x))^2
  return(v)
}


x = matrix(data = h(runif(200*n)), ncol = 200)
estint = apply(x, 2, cumsum) / (1:n)
y = apply(estint, 1, quantile, c(0.025, 0.975))

plot(estint[,1], type = "l", col = 2, ylim = c(0.8, 1.2), xlim = c(0, 10000))
polygon(c(1:n, n:1), c(y[1,], rev(y[2,])), col = "gray90")
lines(estint[,1], type = "l", col = 2, ylim = c(0.8, 1.2))




boot = matrix(data = sample(x[,1], 200*n, replace = TRUE), nrow = n, ncol = 200)
bootit = apply(boot, 2, cumsum) / 1:n
bootU = apply(bootit, 1, quantile, 0.975)
bootL = apply(bootit, 1, quantile, 0.025)


plot(estint[,1], type = "l", col = 2, ylim = c(0.8, 1.2), xlim = c(0, 10000))
lines(bootU, col = 4)
lines(bootL, col = 4)

z = runif(n)
lines(cumsum(h(z))/1:n, col = 5)



plot(x = 1:n, y = theta_est, type = "l", lty=2, ylim = c(0.7, 1.2), xlim = c(0, 10000))
lines(x = 1:n, y = parallelCI[,1], col = "blue")
lines(x = 1:n, y = parallelCI[,3], col = "blue")
lines(x = 1:n, y = cumsum(h(z))/1:n, col = "red", lty = 3)
