#include <RcppArmadillo.h>
#include <cmath> /* pow */
#include <algorithm> /* random_shuffle*/
#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 

arma::vec permDiscoveries(NumericVector ix, NumericVector cv, NumericVector praw) {
  
  int m = ix.length();
  arma::vec one(m);
  one.fill(1);
  //select p-values inside S
  NumericVector PV = praw[ix -1];
  //sort p-values
  std::sort(PV.begin(), PV.end());
  
  //Calculate the number of p-values below cv_1, cv_2, etc as follows:
  int j = 0;
  arma::vec n(m);
  n.fill(0);
  
  for (int i=0; i<m; i++){
    if (PV[i] > cv[j]){
      j = j+1;
    } else{
      n[i] = n[i] + 1;
    }
    if(j >m){
      break;
    }
  }
  
  arma::vec summ(m);
  summ.fill(sum(n));
  arma::vec v(m);
  v = arma::linspace(1, m, m);
  
  arma::vec d(1);
  
  d = max(one - v + summ);
  return (d);
}

/*** R
# library(pARI)
# m <- 1000
# n <- 20
# B <- 100
# delta <- 1
# alpha <- 0.05
# X <- simulateData(0.8,m,n,0,power = 0.8,rho = 0)
# PV <- signTest(X = X,B = 100,seed = rpois(1,1000))
# X<- cbind(PV$pv,PV$pv_H0)
# lambda <- lambdaOpt(X, alpha = alpha, delta = 0, family = "beta")
# cvO <- criticalVector(pvalues = X, lambda = lambda, family = "beta", alpha = alpha,delta = 0)
# ix <- sample.int(m,500)
# max(permDiscoveries1(ix = ix, cv = cvO, praw = PV$pv))
# dI(ix = ix, cv = cvO, pvalues = PV$pv)
*/

// #include <RcppArmadillo.h>
// #include <cmath> /* pow */
// #include <algorithm> /* random_shuffle*/
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::plugins(cpp11)]] 
// // [[Rcpp::export]] 
// 
// int permDiscoveries(NumericVector ix, NumericVector cv, NumericVector praw) {
//   int m = ix.length();
//   NumericVector u (m);
//   for(int i = 1; i<=m; i++){
//     NumericVector PV = praw[ix -1];
//     u[i-1] = 1 - i + sum(PV <= cv[i-1]);
//   }
//   int d = max(u);
//   return (d);
// }

