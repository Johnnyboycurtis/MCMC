library(pixmap)
img <- read.pnm("./img/thorningS.ppm")
plot(img)


img.mat <- img@grey
dim(img.mat)

sdval   <- sd(as.numeric(img.mat))
muvec   <- seq(0.1,0.9,.1)

## ML classification of each pixel
# Homegrown function for maximum likelihood classification
mlclass <- function(img.mat, mean, sd){
  cls <- img.mat
  for (ii in 1:nrow(cls)){
      for (jj in 1:ncol(cls)){
          cls[ii, jj] <- which.max(dnorm(img.mat[ii, jj],
                                          mean, sd, log=T))
      }
  }
  fit   <- cls; fit[] <- mean[cls]
  list(label=cls, fit=fit, img=pixmapGrey(fit))
}


mlc <- mlclass(img.mat, muvec, sdval)

names(mlc)

par(mfrow=c(1,2))
plot(img)
plot(mlc$img)

source("IsingR.R")
sourceCpp("IsingC.cpp")

init <- mlc$label

## R implementation
aa0<-isingR0(img.mat, init, muvec, sdval, beta=4)
plot(aa0$img)

## R implementation slightly faster
aa1<-isingR1(img.mat, init, muvec, sdval, beta=4)
plot(aa1$img)

## C++ implementation
aa2<-isingR2(img.mat, init, muvec, sdval, beta=4)
plot(aa2$img)

## C++ implementation
aa2<-isingR2(img.mat, init, muvec, sdval, beta=4, n_iter=100)
plot(aa2$img)

## Try with random starting image
initr <- sample(seq_along(muvec), length(img.mat), replace=T)
dim(initr) <- dim(img.mat)
aa2<-isingR2(img.mat, initr, muvec, sdval, beta=2, n_iter=500)
plot(aa2$img)

## Try with uniform starting image
initu <- rep(4, length(img.mat))
dim(initu) <- dim(img.mat)
aa2<-isingR2(img.mat, initu, muvec, sdval, beta=2, n_iter=1500)
plot(aa2$img)


## Sampling from the posterior:
isingR3 <- function(img, state, mean, sd, beta=1, n_iter=100, thin=10){
  outiter <- round( n_iter / thin )
  out <- array(NA, c(dim(img), outiter))
  for (k in 1:outiter){
      t0 <- proc.time()
      state <- ising_C(img, state, mean, sd, beta, thin)
      cat(sprintf("time elapsed: %f\n", (proc.time()-t0)[3]))
      out[,,k] <- state
  }
  out
}

aa3<-isingR3(img.mat, init, muvec, sdval, beta=4, n_iter=100)
table(aa3[1,1,])

## Very good question: How to summarize the sequence of images?
## One option is to find the most frequent label of each pixel

i <- 2
j <- 2
aa3[i,j,]
table( aa3[i,j,] )

map <- rep(NA, length(img.mat))
dim(map) <- dim(img.mat)
for (i in 1:nrow(map)){
    for (j in 1:ncol(map)){
        map[i,j] <- as.numeric( names(which.max(table(aa3[i,j,]))))
    }
}

fit <- muvec[map]
dim(fit) <- dim(img.mat)
img <- pixmapGrey(fit)


## Do SVD on image

S<-svd(img.mat)
d <- S$d
v <- S$v
u <- S$u

d2<- d
d2[-(1:150)] <- 0
R<-u%*%(d2*t(v))
par(mfrow=c(1,2))
plot(pixmapGrey(img.mat))
plot(pixmapGrey(R))




