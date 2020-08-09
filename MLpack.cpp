// [[Rcpp::depends(RcppMLPACK)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppMLPACK.h>
#include <mlpack/methods/linear_regression/linear_regression.hpp>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
List kmeans_cpp(const arma::mat &data,
                const int &clusters) {
  // ...
  arma::Col<size_t> assignments;
  mlpack::kmeans::KMeans<> k;
  k.Cluster(data, clusters, assignments);
  
  return List::create(
    _["clusters"] = clusters,
    _["result"] = assignments
  );
}

// [[Rcpp::export]]
Rcpp::List logisticRegression(const arma::mat& train, 
                              const arma::vec& labels,
                              const Rcpp::Nullable<Rcpp::NumericMatrix>& test = R_NilValue) {
  // MLPACK wants Row<size_t> which is an unsigned representation that R does not have
  arma::vec resultsur;
  // TODO: check that all values are non-negative
  // labelsur = arma::conv_to<arma::Col<size_t>>::from(labels);
  // Initialize with the default arguments. TODO: support more arguments>
  mlpack::regression::LogisticRegression<> lrc(train, labels);
  arma::rowvec parameters = lrc.Parameters();
  Rcpp::List return_val;
  if (test.isNotNull()) {
    arma::mat test2 = Rcpp::as<arma::mat>(test);
    lrc.Predict(test2, resultsur);
    arma::vec results = arma::conv_to<arma::vec>::from(resultsur);
    return_val = Rcpp::List::create(
      Rcpp::Named("parameters") = parameters,
      Rcpp::Named("results") = results
    );
  } 
  else {
    return_val = Rcpp::List::create(Rcpp::Named("parameters") = parameters);
  }
  return return_val;
}