#include <RcppArmadillo.h>
#include <math.h>
#include <cmath> /* pow */
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 


NumericVector lambdaCalibrate(arma::mat X, arma::vec alpha, int delta, std::string family) {
  int m = X.n_rows;
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
    if(family == "finner"){

     lambda = ((mV - 1)%Y.rows(delta, m-1).col(bb))/((1-Y.rows(delta, m-1).col(bb))%((idV-deltaV)*alpha));
    }
    if(family == "higher.criticism"){

     lambda = (sqrt(m)*((idV/m) - Y.col(bb)))/(sqrt(Y.col(bb)%(1-Y.col(bb))));
      
    }
    if(family == "beta"){

      for (int i=1; i<m; i++) {
      double q = arma::conv_to<double>::from(Y.col(bb).row(i));
      double shape1 = i;
      double shape2 = m+1-i;
      long double beta = 1-exp(R::pbeta(q=q,shape1,shape2,0,1));
      if (beta == 0){
        beta = 1;
      }
      lambda[i-1] = beta;
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
library(ARIpermutation)
#m <- 1000
#n <- 20
#B <- 100
#delta <- 1
#alpha <- 0.05
#X <- simulateData(0.9,m,n,0,power = 0.8)
#PV <- signTest(X = t(X),B = 100)
#X<- cbind(PV$pv,PV$pv_H0)
#Test <- lambdaCalibrate(X, alpha = alpha, delta = 0, family = "beta")
#Test
#lambdaOpt(t(X), family = "beta", alpha = 0.05, delta =0)
*/

