#include <RcppArmadillo.h>
#include <math.h>
#include <cmath> /* pow */
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 


NumericVector lambdaCalibrate(arma::mat X, arma::vec alpha, double delta, std::string family) {
  double m = X.n_rows;
  int B = X.n_cols;
  NumericVector T(B);
  NumericVector lambda(m-delta);
  NumericVector idx(1);
  arma::vec idV(m-delta);
  std::iota(idV.begin(), idV.end(), 1 + delta);
  arma::vec deltaV(m-delta);
  deltaV.fill(delta);
  arma::vec mV(m-delta);
  mV.fill(m);
  
  //sort columns ascending order
  arma::mat Y(m, B, arma::fill::zeros);
  for (int rr=0; rr<B; rr++) {
    arma::colvec x = X.col(rr); 
    std::sort(x.begin(), x.end());
    Y.col(rr) = x;
  }
  
//  std::iota(idV.begin(), idV.end(), 1);
  
  //compute lambda for each permutations
  for (int bb=0; bb<B; bb++) {
    if(family == "simes"){
 //    int minc = delta;
  //   int maxc = m;
     lambda = (mV-deltaV)%(Y.rows(delta, m-1).col(bb))/((idV-deltaV)*alpha);

    }
    if(family == "aorc"){

     lambda = ((mV - 1)%Y.rows(delta, m-1).col(bb))/((1-Y.rows(delta, m-1).col(bb))%((idV-deltaV)*alpha));
    }
    if(family == "higher.criticism"){

     lambda = (sqrt(m)*((idV/m) - Y.col(bb)))/(sqrt(Y.col(bb)%(1-Y.col(bb))));
      
    }
    if(family == "beta"){

      for (int i=0; i<m; i++) {
      long double q = arma::conv_to<double>::from(Y.col(bb).row(i));
      double shape1 = i+1;
      double shape2 = m-i;
     // long double beta = 1-exp(R::pbeta(q=q,shape1,shape2,1,1));
      lambda[i] = R::pbeta(q = q, shape1, shape2, 1, 0);
      }
    }
//take minimum over hypotheses for each permutations
    T[bb] = min(lambda);
    
  }
  std::sort(T.begin(), T.end());
  idx = floor(alpha*B);
//  lambdaE = T(floor(alpha*B)+1);
  return (T[idx]);
}


/*** R
# library(pARI)
# m <- 1000
# n <- 20
# B <- 100
# delta <- 1
# alpha <- 0.05
# X <- simulateData(0.9,m,n,0,power = 0.8,set.seed = rpois(1,1000))
# PV <- signTest(X = X,B = 100,seed = rpois(1,1000))
# X<- cbind(PV$pv,PV$pv_H0)
# Test <- lambdaCalibrate(t(X), alpha = alpha, delta = 0, family = "simes")
# Test
# lambdaOpt1(X, family = "simes", alpha = alpha, delta =0)
*/

