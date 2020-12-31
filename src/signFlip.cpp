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
  arma::vec eps, M, Var, X2, M2, id;
  
  //X = X / sqrt(n);    // scaling
  
  int bb;
  for (bb=0; bb<B; bb++) {
    eps = Rcpp::rbinom(n, 1, 0.5)*2 - 1;  
    id = Rcpp::rbinom(n, 1, 1)*2 - 1;   
    M = (X * eps) / n; 
    X2 = (pow(X, 2) * id)/n; 
    M2 = pow(M,2)/pow(n,2); 
    Var = (X2 - M2)*(n/(n-1));
    T.col(bb) = M / sqrt(Var/n);
    
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

