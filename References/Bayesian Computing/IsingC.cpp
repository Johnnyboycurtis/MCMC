//#####################################################################
//## Metropolis sampler for Ising model                              ##
//##                                                                 ##
//## Author: Søren Højsgaard, May 2014                               ##
//#####################################################################

#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
IntegerMatrix ising2_C(NumericMatrix imgmat, IntegerMatrix init, NumericVector mean,
             double sd, double beta, int n_iter){

  int nr=imgmat.rows(), nc=imgmat.cols(), nn=imgmat.rows()*imgmat.cols();
  int ii, jj, uu, cnb, pnb, cij, pij, nclass = mean.length();
  double hr, alpha;

  IntegerVector labels    = Range(1, nclass);
  IntegerMatrix state     = IntegerMatrix(nr, nc, init.begin());
  NumericMatrix curr_logL = NumericMatrix(nr, nc);
  NumericMatrix prop_logL = NumericMatrix(nr, nc);

  for (ii=0; ii<nn; ii++){
    curr_logL[ii] = -0.5*pow((imgmat[ii]-mean[state[ii]-1])/sd, 2);
  }

  for (uu=0; uu<n_iter; uu++){
    IntegerVector NN   = Rcpp::RcppArmadillo::sample(labels, nn, true);
    IntegerMatrix pmat = IntegerMatrix(nr, nc, NN.begin());
    arma::mat unifmat  = randu<arma::mat>(nr, nc);

    for (ii=0; ii<nn; ii++){
      prop_logL[ii] = -0.5*pow((imgmat[ii]-mean[pmat[ii]-1])/sd, 2);
    }

    for (ii=1; ii<(nr-1); ii++){
      for (jj=1; jj<(nc-1); jj++){
        cij= state(ii,jj);
        pij= pmat(ii,jj);
        cnb = (cij==state(ii-1,jj)) +  (cij==state(ii+1,jj)) +
          (cij==state(ii,jj-1)) +  (cij==state(ii,jj+1));
        pnb = (pij==state(ii-1,jj)) +  (pij==state(ii+1,jj)) +
          (pij==state(ii,jj-1)) +  (pij==state(ii,jj+1));
        hr  = exp( prop_logL(ii,jj)-curr_logL(ii,jj)+beta*(pnb-cnb) );
        if (hr < 1) alpha=hr; else alpha=1;
        if (alpha > unifmat(ii,jj)) {
          state(ii, jj)     = pmat(ii, jj);
          curr_logL(ii, jj) = prop_logL(ii, jj);
        }
      }
    }
  }
  return state;
}


/*** R

isingR2 <- function(img, state, mean, sd, beta=1, n_iter=10){
  t0 <- proc.time()
  state <- ising2_C(img, state, mean, sd, beta, n_iter)
  cat(sprintf("time elapsed: %f\n", (proc.time()-t0)[3]))
  fit  	   <- muvec[state]
  dim(fit) <- dim(state)
  list(label=state, fit=fit, img=pixmapGrey(fit))
}

***/

