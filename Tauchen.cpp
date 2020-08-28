#include <RcppArmadillo.h>
#include <Tauchen.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List tauchen_test(const double rho,
                  const double sigmu_u,
                  const int m=3,
                  const int n=7) {
  // ...
  Tauchen tc(rho, sigmu_u, 0.0, m, n);
  return List::create(
    _["grid"] = tc.state_values,
    _["pi"] = tc.P,
    _["station"] = tc.stationary_distributions
  );
}