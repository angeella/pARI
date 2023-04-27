#include <RcppArmadillo.h>
#include <math.h>
#include <cmath> /* pow */
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 


NumericVector lambdaCalibrate(arma::mat X, arma::vec alpha, double delta, std::string family, double m) {
  int B = X.n_cols;
  int mm = X.n_rows;
  NumericVector T(B);
  NumericVector lambda(mm-delta);
  NumericVector idx(1);
  arma::vec idV(mm-delta);
  std::iota(idV.begin(), idV.end(), 1 + delta);
  arma::vec deltaV(mm-delta);
  deltaV.fill(delta);
  arma::vec mV(mm-delta);
  mV.fill(m);
  
  //sort columns ascending order
  arma::mat Y(mm, B, arma::fill::zeros);
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
     lambda = (mV-deltaV)%(Y.rows(delta, mm-1).col(bb))/((idV-deltaV)*alpha);

    }
    if(family == "aorc"){

     lambda = ((mV - 1)%Y.rows(delta, mm-1).col(bb))/((1-Y.rows(delta, mm-1).col(bb))%((idV-deltaV)*alpha));
    }
    if(family == "higher.criticism"){

     lambda = (sqrt(mV)%((idV/mV) - Y.col(bb)))/(sqrt(Y.col(bb)%(1-Y.col(bb))));
      
    }
    if(family == "beta"){

      for (int i=0; i<mm; i++) {
      long double q = arma::conv_to<double>::from(Y.col(bb).row(i));
      double shape1 = i+1;
      double shape2 = m-i;
     // long double beta = 1-exp(R::pbeta(q=q,shape1,shape2,1,1));
      lambda[i] = R::pbeta(q, shape1, shape2, 1, 0);
      }
    }
    if(family == "power"){
      lambda = log(Y.col(bb))/log(idV/mV)*alpha;
    }
//take minimum over hypotheses for each permutations
    T[bb] = min(lambda);
    
  }
  std::sort(T.begin(), T.end());
  idx = floor(alpha*B);
//  lambdaE = T(floor(alpha*B)+1);
  
  return (T[idx]);
//  return(quantile(T, alpha, 1));
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
# Test <- lambdaCalibrate(X, alpha = alpha, delta = 0, family = "simes")
# Test
# lambdaOptR(t(X), family = "simes", alpha = alpha, delta =0)
*/

