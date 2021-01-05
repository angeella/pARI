#include <RcppArmadillo.h>
#include <cmath> /* pow */
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 
arma::mat signFlip(arma::mat X, double B) {
  double m = X.n_rows;
  double n = X.n_cols;
  
  
  arma::mat T(m, B, arma::fill::zeros);
  arma::vec eps, Tb, Tb0, Tb1, Tb2, eps1;
  
  //X = X / sqrt(n);    // scaling
  
  int bb;
  for (bb=0; bb<B; bb++) {
    eps = Rcpp::rbinom(n, 1, 0.5)*2 - 1;  // signs 
    eps1 = Rcpp::rbinom(n, 1, 1)*2 - 1;  // identity    
    Tb = X * eps;
    Tb1 = Tb / n; //mean
    Tb0 = (pow(X, 2) * eps1)/n; //E(x^2)
    Tb2 = pow(Tb1,2); //E(X)^2
    Tb =  (Tb0-Tb2)*(n/(n-1)); //sample var
    T.col(bb) = Tb1/sqrt(Tb/n);
  }
  return (T);
}


/*** R
#m <- 100
#n <- 100
#B <- 1000
#set.seed(123)
#X <- matrix(rnorm(m*n), ncol=n)
#set.seed(123)
#T <- signFlip(X, B)
#str(T)
*/

