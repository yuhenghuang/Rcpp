#include <algorithm>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
void stl_sort_inplace(arma::vec &x) {
  std::sort(x.begin(), x.end());
}


// [[Rcpp::export]]
arma::vec stl_sort(arma::vec &x) {
  // copy x = y
  arma::vec y = x;
  std::sort(y.begin(), y.end());
  return y;
}
