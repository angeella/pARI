#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Sorting of the rows of a matrix by ascending order
//
// @param X A numeric matrix.
// @return A numeric matrix whose rows are sorted by ascending order.
//
// @export
// [[Rcpp::export]]
arma::mat rowSortC(arma::mat X) {
  // sorting by row in ascending order
  int k = X.n_rows;
  int B = X.n_cols;
  arma::mat Y(k, B, arma::fill::zeros);
  for (int rr=0; rr<k; rr++) {
    arma::rowvec x = X.row(rr);  // ascending order
    std::sort(x.begin(), x.end());
    Y.row(rr) = x;
  }
  return Y;
}

/*** R
A <- matrix(rnorm(15), 5, 3);
B <- rowSortC(A);
*/

