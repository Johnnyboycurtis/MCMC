#####################################################################
## Metropolis sampler for Ising model                              ##
##                                                                 ##
## Author: Søren Højsgaard, May 2014                               ##
#####################################################################

library(pixmap); library(compiler); enableJIT(3)

isingR0 <- function(img, state, mean, sd, beta=1, n_iter=10){
  t0 <- proc.time()
  state <- ising_R0(img, state, mean, sd, beta, n_iter)
  cat(sprintf("time elapsed: %f\n", (proc.time()-t0)[3]))
  fit  	   <- muvec[state]
  dim(fit) <- dim(state)
  list(label=state, fit=fit, img=pixmapGrey(fit))
}

ising_R0 <- function(img, state, mean, sd, beta=1, n_iter=10){
  curr.logL      <- dnorm(img, mean=mean[state], sd=sd, log=TRUE)
  dim(curr.logL) <- dim( state )
  label <- seq_along(mean)
  for (uu in 1:n_iter){
      state.old      <- state
      prop           <- sample(label, length(state), replace=TRUE)
      dim(prop)      <- dim( state )
      prop.logL      <- dnorm(img, mean=mean[prop], sd=sd, log=TRUE)
      dim(prop.logL) <- dim( state )
      unif           <- runif(length(state))
      dim(unif)      <- dim( state )
      for (ii in 2:(nrow( state )-1)){      # Visit each site (i,j)
          for (jj in 2:(ncol( state )-1)){  # ignoring boundaries
              idxmat <- matrix(c(ii-1,ii+1,ii,ii,  jj,jj,jj-1,jj+1), nrow=4)
              sc.ij  <- state[ii,jj] # current state
              sp.ij  <- prop[ii,jj]  # proposed state
              nb     <- state[idxmat]
              nc.ij  <- sum(sc.ij == nb)
              np.ij  <- sum(sp.ij == nb)
              hr     <- exp( (prop.logL[ii,jj]-curr.logL[ii,jj] + beta*(np.ij-nc.ij)) )
              alpha  <- min( 1, hr )
              if (alpha > unif[ii,jj]) {
                  state[ ii, jj ]     <- sp.ij
                  curr.logL[ ii, jj ] <- prop.logL[ ii, jj ]
              }
          }
      }
      new.states <- sum((state-state.old)!=0)
      cat(sprintf("iteration %3i new states: %6d pct. new states: %5.2f\n",
                  uu, new.states, 100*new.states/length(state)))
  }
  state
}


isingR1 <- function(img, state, mean, sd, beta=1, n_iter=10){
  t0 <- proc.time()
  state <- ising_R1(img, state, mean, sd, beta, n_iter)
  cat(sprintf("time elapsed: %f\n", (proc.time()-t0)[3]))
  fit  	   <- muvec[state]
  dim(fit) <- dim(state)
  list(label=state, fit=fit, img=pixmapGrey(fit))
}


ising_R1 <- function(img, state, mean, sd, beta=1, n_iter=10){
  curr.logL      <- dnorm(img, mean=mean[state], sd=sd, log=TRUE)
  dim(curr.logL) <- dim( state )
  for (uu in 1:n_iter){
    prop           <- sample(seq_along(mean), length(state), replace=TRUE)
    dim(prop)      <- dim( state )
    prop.logL      <- dnorm(img, mean=mean[prop], sd=sd, log=TRUE)
    dim(prop.logL) <- dim( state )
    unif           <- runif(length(state))
    dim(unif)      <- dim( state )
    idxmat         <- matrix(c(0,1, 2,1, 1,0, 1,2), nrow=4, byrow=TRUE)
    state.old      <- state
    for (ii in 2:(nrow( state )-1)){  # Visit each site (i,j)
      idxmat[, 1] <- idxmat[,1] + c(1,1,1,1)
      idxmat[, 2] <- c(1,1,0,2)
      for (jj in 2:(ncol( state )-1)){    # ignoring boundaries
        sc.ij  <- state[ii,jj] # current state
        sp.ij  <- prop[ii,jj]  # proposed state
        idxmat <- idxmat + c(0,0,0,0,1,1,1,1)
        nb     <- state[idxmat]
        nc.ij  <- sum(sc.ij == nb)
        np.ij  <- sum(sp.ij == nb)
        hr     <- exp( (prop.logL[ii,jj]-curr.logL[ii,jj] + beta*(np.ij-nc.ij)) )
        alpha  <- min( 1, hr)
        if(alpha > unif[ii,jj]) {
          state[ ii, jj ]     <- sp.ij
          curr.logL[ ii, jj ] <- prop.logL[ ii, jj ]
        }
      }
    }
    new.states <- sum((state-state.old)!=0)
    cat(sprintf("iteration %3i new states: %6d pct. new states: %5.2f\n",
                uu, new.states, 100*new.states/length(state)))
  }
  state
}

