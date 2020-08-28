#include <RcppArmadillo.h>
#include <CubicSpline.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

// [[Rcpp::export]]
arma::vec cubic_spline_test(arma::vec &X,
                            arma::vec &Y,
                            arma::vec &t) {
  // ...
  cubic_spline cs(X.begin(), Y.begin(), X.n_elem);

  arma::vec res = cs(t);

  return res;
}