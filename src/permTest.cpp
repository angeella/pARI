#include <RcppArmadillo.h>
#include <cmath> /* pow */
#include <algorithm> /* random_shuffle*/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 
arma::mat permT(arma::mat X, double B, arma::vec label) {
  int m = X.n_rows;

  double n1, n2;
  arma::vec permLabel, M1, M2, Tb01, Tb02, Tb1, Tb2, Tb21, Tb22, pV, I1, I2;
  std::vector<int>::iterator id;

  n1 = std::count(label.begin(), label.end(), 1);
  n2 = std::count(label.begin(), label.end(), 2);
  arma::mat T(m, B, arma::fill::zeros);
  arma::mat X1(m, n1, arma::fill::zeros);
  arma::mat X2(m, n2, arma::fill::zeros);
  I1 = Rcpp::rbinom(n1, 1, 1)*2 - 1;
  I2 = Rcpp::rbinom(n2, 1, 1)*2 - 1;
  
  int bb;
  for (bb=0; bb<B; bb++) {
    
    std::random_device rd;
    std::mt19937 g(rd());
    
    std::shuffle(label.begin(), label.end(), g);
    
    X1 = X.cols(find(label == 1));
    X2 = X.cols(find(label == 2)); 
    M1 = (X1 * I1)/n1;
    M2 = (X2 * I2)/n2;
    Tb01 = (pow(X1, 2) * I1)/n1; //E(x^2)
    Tb21 = pow(M1,2)/pow(n1,2); //E(X)^2
    Tb1 =  (Tb01-Tb21)*(n1/(n1-1)); //sample var
    Tb02 = (pow(X2, 2) * I2)/n2; //E(x^2)
    Tb22 = pow(M2,2)/pow(n2,2); //E(X)^2
    Tb2 =  (Tb02-Tb22)*(n2/(n2-1)); //sample var
 //   pV =  ((n1 - 1)* Tb1 + (n2 - 1)* Tb2)/ (n1 + n2 - 2);
    pV = Tb1/n1 + Tb2/n2;
    T.col(bb) = (M1 - M2)/sqrt(pV * (1/n1 + 1/n2));
   // T.col(bb) = (M1 - M2)/sqrt(pV);
   // T.col(bb) = (M1 - M2);
   // T.col(bb) = pV * (1/n1 + 1/n2);
  }
  return (T);
}


/*** R
#library(multtest)
#data(golub)
#X = golub
#B = 100
#label = factor(golub.cl)
#Test<- permTest(X , B, label)
#str(Test)
*/