#include <RcppArmadillo.h>
#include <cmath> /* pow */
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 
arma::mat signFlip(arma::mat X, double B) {
  int m = X.n_rows;
  int n = X.n_cols;
  
  arma::mat T(m, B, arma::fill::zeros);
  
  arma::vec eps, M, id, Tb2, V, T0, Tb;
  
  int bb;
  for (bb=0; bb<B; bb++) {
    eps = Rcpp::rbinom(n, 1, 0.5)*2 - 1;  // signs
    id = Rcpp::rbinom(n, 1, 1)*2 - 1; //identity
    Tb = X * eps;
    M = Tb /n; //mean
    //    Tb2 = pow (Tb,2)/pow (n,2); //mean power 2
    //    Tbp = (pow (X,2) * id)/n; //power 2 mean
    Tb2 = pow(X, 2) * id; 
    //Tb0 = Tb0/n;
    V = (Tb2 - n * pow(M,2))/(n-1);
    //    Ts = Tb/ pow (((Tbp-Tb2)*n)/(n-1), 0.5); //final test
    T0 = M/ pow((V)/n,0.5);
    T.col(bb) = T0;
    
  }
  return (T);
}


/*** R
#m <- 100
#n <- 10
#B <- 200
#X <- matrix(rnorm(m*n), ncol=n)
#set.seed(123)
#T <- signFlip(X, B)
#str(T)
*/





arma::mat testByPermutation(arma::mat X, NumericVector cls, double B) {
  int m = X.n_rows;
  //    int n = X.n_cols;
  arma::mat T(m, B, arma::fill::zeros);
  
  // arma::vec eps;
  // arma::mat Tb;
  //    RNGScope scope;
  int bb, ii;
  for (bb=0; bb<B; bb++) {
    for (ii=0; ii<m; ii++) {
      // TODO!
    }
  }
  return(T);
}