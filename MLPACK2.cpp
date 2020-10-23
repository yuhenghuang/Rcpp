#include <RcppMLPACK.h>
#include <mlpack/methods/linear_regression/linear_regression.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>
#include <mlpack/core/tree/cover_tree.hpp>
#include <mlpack/core/metrics/lmetric.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
// [[Rcpp::depends(RcppMLPACK)]]


using namespace Rcpp;


// [[Rcpp::export]]
arma::vec linear_regression(const arma::mat& X,
                            const arma::vec& Y,
                            const double lambda = 0.0,
                            const bool intercept = true) {
  // ...

  Rcout << "reach here" << std::endl;

  Rcout << X.n_rows << ", " << X.n_cols << std::endl;
  Rcout << Y.n_rows << ", " << Y.n_cols << std::endl;
  mlpack::regression::LinearRegression lr(X, Y.t(), lambda, intercept);


  arma::rowvec Ybar(Y.n_elem);
  lr.Predict(X, Ybar);
  return Ybar.t();
}


// [[Rcpp::export]]
List kmeans(const arma::mat& data,
            const int clusters) {
  // ...

  Rcout << data.n_rows << ", " << data.n_cols << std::endl;
  arma::Row<size_t> assignments;
  arma::mat centroids(data.n_rows, clusters);
  mlpack::kmeans::KMeans<> k;
  k.Cluster(data, clusters, assignments, centroids);

  Rcout << centroids.n_rows << ", " << centroids.n_cols << std::endl;

  return List::create(
    _["clusters"] = clusters,
    _["assignments"] = assignments,
    _["centroids"] = centroids
  );
}


// [[Rcpp::export]]
List logistic_regression(const arma::mat& train,
                         const arma::irowvec& labels,
                         const Nullable<NumericMatrix>& test = R_NilValue) {
  // ...
  arma::Row<size_t> labelsur, resultsur;

  labelsur = arma::conv_to<arma::Row<size_t>>::from(labels);

  mlpack::regression::LogisticRegression<> lrc(train, labelsur);

  arma::rowvec parameters = lrc.Parameters();

  List return_val;

  if (test.isNotNull()) {
    arma::mat test2 = as<arma::mat>(test);
    lrc.Classify(test2, resultsur);
    arma::vec results = arma::conv_to<arma::vec>::from(resultsur);
    return_val = List::create(_["parameters"] = parameters,
                              _["results"] = results);
  }
  else {
    return_val = List::create(_["parameters"] = parameters);
  }

  return return_val;
}


// [[Rcpp::export]]
List naive_bayes_classifier(const arma::mat& train,
                            const arma::irowvec& labels,
                            const int classes,
                            const Nullable<NumericMatrix>& test = R_NilValue) {
  // ...

  arma::Row<size_t> labelsur, resultsur;

  labelsur = arma::conv_to<arma::Row<size_t>>::from(labels);

  mlpack::naive_bayes::NaiveBayesClassifier<> nbc(train, labelsur, classes);

  List return_val;
  if (test.isNotNull()) {
    arma::mat test2 = as<arma::mat>(test);
    nbc.Classify(test2, resultsur);

    arma::irowvec results = arma::conv_to<arma::irowvec>::from(resultsur);

    return_val = List::create(_["means"] = nbc.Means(),
                              _["vars"] = nbc.Variances(),
                              _["probs"] = nbc.Probabilities(),
                              _["classifications"] = results);
  }
  else {
    return_val = List::create(_["means"] = nbc.Means(),
                              _["vars"] = nbc.Variances(),
                              _["probs"] = nbc.Probabilities());
  }

  return return_val;
}


// [[Rcpp::export]]
List nn_ml(const arma::mat &data,
           const int k) {
  // ...
  mlpack::tree::StandardCoverTree<mlpack::metric::EuclideanDistance, 
                                  mlpack::neighbor::NeighborSearchStat<mlpack::neighbor::NearestNeighborSort>,
                                  arma::mat> tree(data);

  mlpack::neighbor::NeighborSearch<mlpack::neighbor::NearestNeighborSort,
                                   mlpack::metric::LMetric<2>,
                                   arma::mat,
                                   mlpack::tree::StandardCoverTree> coverTreeSearch(
                                     std::move(tree),
                                     mlpack::neighbor::SINGLE_TREE_MODE, 
                                     0.05);

  arma::Mat<size_t> neighborsCoverTree;
  arma::mat distanceCoverTree;
  coverTreeSearch.Search(data, k, neighborsCoverTree, distanceCoverTree);
  
  return List::create(
    _["clusters"] = neighborsCoverTree,
    _["result"] = distanceCoverTree
  );
}
