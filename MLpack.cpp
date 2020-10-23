#include <RcppMLPACK.h>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/core/tree/cover_tree.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>

// [[Rcpp::depends(RcppMLPACK)]]
using namespace Rcpp;

// [[Rcpp::export]]
List kmeans_ml(const arma::mat &data,
							 const int &clusters) {
	// ...
	arma::Row<size_t> assignments;
	mlpack::kmeans::KMeans<> k;
	k.Cluster(data, clusters, assignments);
  
	return List::create(
		_["clusters"] = clusters,
		_["result"] = assignments
	);
}

// [[Rcpp::export]]
List logistic_regression(const arma::mat &train, 
                         const arma::irowvec &labels,
                         const Nullable<NumericMatrix> &test = R_NilValue) {
  // MLPACK wants Row<size_t> which is an unsigned representation that R does not have
  arma::Row<size_t> labelsur, resultsur;
  // TODO: check that all values are non-negative
  labelsur = arma::conv_to<arma::Row<size_t>>::from(labels);

  // Initialize with the default arguments. TODO: support more arguments>
  mlpack::regression::LogisticRegression<> lrc(train, labelsur);
  arma::rowvec parameters = lrc.Parameters();
  List return_val;
  if (test.isNotNull()) {
    arma::mat test2 = Rcpp::as<arma::mat>(test);
    lrc.Classify(test2, resultsur);
    arma::vec results = arma::conv_to<arma::vec>::from(resultsur);
    return_val = Rcpp::List::create(
      Rcpp::Named("parameters") = parameters,
      Rcpp::Named("results") = results
    );
  } 
  else {
    return_val = List::create(Rcpp::Named("parameters") = parameters);
  }
  return return_val;
}


using namespace mlpack::neighbor;
using namespace mlpack::tree;
using namespace mlpack::metric;
// [[Rcpp::export]]
List nn_ml(const arma::mat &data,
           const int k) {
  // ...
  CoverTree<LMetric<2, true>,
            NeighborSearchStat<NearestNeighborSort>,
            arma::mat> tree(data);

  NeighborSearch<NearestNeighborSort, 
                 LMetric<2, true>,
                 arma::mat,
                 mlpack::tree::KDTree<LMetric<2, true>,
                                      NeighborSearchStat<NearestNeighborSort>,
                                      arma::mat>> cover_tree_search(tree);
                                  
  arma::Mat<size_t> cover_tree_neighbors;
  arma::mat cover_tree_distance;
  cover_tree_search.Search(k, cover_tree_neighbors, cover_tree_distance);

  return List::create(
    _["clusters"] = cover_tree_neighbors,
    _["result"] = cover_tree_distance
  );
}
