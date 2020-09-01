#define ARMA_USE_SUPERLU 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
void superLU_demo() {
  arma::sp_mat A = arma::sprandu<arma::sp_mat>(1000, 1000, 0.1);

  arma::vec b = arma::randu<arma::vec>(1000);
  arma::mat B = arma::randu<arma::mat>(1000, 5);

  arma::vec x = arma::spsolve(A, b);
  arma::mat X = arma::spsolve(A, B);

  bool status = arma::spsolve(x, A, b);
  if (!status)
    Rcout << "no solution" << std::endl;

  arma::spsolve(x, A, b, "lapack");
  arma::spsolve(x, A, b, "superlu");

  Rcout << "Done." << std::endl;
}
