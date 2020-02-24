#include <RcppArmadillo.h>
#include <cmath> /* pow */
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 

// [[Rcpp::export]] 

arma::mat varBySignFlipping(arma::mat X, double B) {
    int m = X.n_rows;
    int n = X.n_cols;

   
    arma::mat T(m, B, arma::fill::zeros);
    arma::vec eps, Tb, Tb0, Tb1, eps1;
    
    //X = X / sqrt(n);    // scaling
    
    int bb;
    for (bb=0; bb<B; bb++) {
        eps = Rcpp::rbinom(n, 1, 0.5)*2 - 1;  // signs 
        eps1 = Rcpp::rbinom(n, 1, 1)*2 - 1;  // identity    
        Tb = X * eps;
        Tb1 = Tb / n;
        Tb0 = pow(X, 2) * eps1;
        //Tb0 = Tb0/n;
        Tb = (Tb0 - n * pow(Tb1,2))/(n-1); //(without Tb0 = Tb0/n)
        //Tb = (Tb0 - pow(Tb1,2));
        T.col(bb) = Tb;
        
    }
    
   return(T);
}