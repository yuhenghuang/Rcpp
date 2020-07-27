#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
void showValue(double x) {
  Rcout << "The value is " << x << std::endl;
}

// [[Rcpp::export]]
void showMatrix(arma::mat X) {
  Rcout << "Armadillo matrix is" << '\n' << X << std::endl;
}

// [[Rcpp::export]]
void callPrint(RObject x) {
  print(x);
}

// [[Rcpp::export]]
void useOperatorOnVector(NumericVector x) { 
  // C print
  Rcout << "Rcpp vector is " << '\n' << x << std::endl;
  // call R internal print
  print(x);
}

// [[Rcpp::export]]
void useOperatorOnMatrix(NumericMatrix x) { 
  Rcout << "Rcpp matrix is " << '\n' << x << std::endl;
  print(x);
}