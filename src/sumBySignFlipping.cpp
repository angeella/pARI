#include <RcppArmadillo.h>
#include <cmath> /* pow */
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 

// [[Rcpp::export]] 


arma::mat meanBySignFlipping(arma::mat X, double B) {
    int m = X.n_rows;
    int n = X.n_cols;

    arma::mat T(m, B, arma::fill::zeros);
    
    arma::vec eps, Tb;
    
    //X = X / sqrt(n);    // scaling
    
    int bb;
    for (bb=0; bb<B; bb++) {
        eps = Rcpp::rbinom(n, 1, 0.5)*2 - 1;  // signs
        Tb = X * eps;
        Tb = Tb /n;
        T.col(bb) = Tb;

    }
   return(T);
}