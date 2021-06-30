#include <RcppArmadillo.h>
#include <cmath> /* pow */
#include <algorithm> /* random_shuffle*/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 

int permDiscoveries(NumericVector ix, NumericVector cv, NumericVector praw) {
  int m = ix.size();
  NumericVector u(m);
  for(int i = 1; i<=m; i++){
    NumericVector PV = praw[ix -1];
    u[i-1] = 1 - i + sum(PV <= cv[i-1]);
  }
  int d = max(u);
  return (d);
}

/*** R
library(pARI)
#m <- 1000
#n <- 20
#B <- 100
#delta <- 1
#alpha <- 0.05
#X <- simulateData(0.9,m,n,0,power = 0.8,set.seed = rpois(1,1000))
#PV <- signTest(X = t(X),B = 100,seed = rpois(1,1000))
#X<- cbind(PV$pv,PV$pv_H0)
#lambda <- lambdaOpt(X, alpha = alpha, delta = 0, family = "beta")
#cvO <- cv(pvalues = X, lambda = lambda, family = "beta", alpha = alpha, ct = c(0,1),delta = 0)
#ix <- sample.int(m,500)
#permDiscoveries(ix = ix, cv = cvO, praw = PV$pv)
#dI(ix = ix, cv = cvO, praw = PV$pv)
*/