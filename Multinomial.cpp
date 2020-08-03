#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

arma::Col<int> rmultinom_1(const uword &size,
                           arma::vec &probs,
                           const uword &N) {
  // ...
  arma::Col<int> out(N);
  rmultinom(size, probs.begin(), N, out.begin());
  return out;
}

// [[Rcpp::export]]
arma::Mat<int> rmultinom_arma(const uword &n,
                              const uword &size,
                              arma::vec &probs) {
  // ...
  uword N = probs.size();
  arma::Mat<int> sim(N, n);
  for (uword i=0; i<n; ++i)
    sim.col(i) = rmultinom_1(size, probs, N);
  return sim;
}

// [[Rcpp::export]]
arma::Mat<int> rmultinom_omp(const uword &n,
                             const uword &size,
                             arma::vec &probs,
                             const int &cores=1) {
  // ...
  omp_set_num_threads(cores);

  uword N = probs.size();
  arma::Mat<int> sim(N, n);

  #pragma omp parallel for
  for (uword i=0; i<n; ++i)
    sim.col(i) = rmultinom_1(size, probs, N);
  return sim;
}

// [[Rcpp::export]]
IntegerVector rpois_arma(const uword &n,
                          const arma::vec &lambda) {
  // ...
  uword size = lambda.size();
  IntegerVector sim(n);
  for (uword i=0; i<n; ++i) 
    sim[i] = R::rpois(lambda[i%size]);
  return sim;
}