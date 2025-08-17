// [[Rcpp::depends(RcppArmadillo)]]
/*
  // [[Rcpp::plugins(cpp11)]]
*/

#include <RcppArmadillo.h>
#include <omp.h>
#include <vector>
#include <algorithm>


using namespace Rcpp;
using arma::uword;

const static int NUM_THREAD = 10;

//' Regression with cluster
//'
//' @param X (n, k) matrix (double), independent variable
//' @param Y (n, 1) column vector (double), dependent variable
//' @param G (n, 1) column vector (unsigned integer), grouping flags
//'
//' @return List of 
//' parameters: (k, 1) arma::vec,
//' robust std.err: (k, 1) arma::vec
//'
// [[Rcpp::export]]
List cluster_regression(const arma::mat &X,
                        const arma::vec &Y,
                        const arma::uvec &G) {
  // ...
  uword n = X.n_rows, k = X.n_cols, c = G.max()+1;
  double df=n-k;

  arma::mat outer = (X.t() * X).i();
  arma::vec beta = outer * X.t() * Y;
  arma::vec e = Y - X * beta;

  arma::ucolvec G_count(c, arma::fill::zeros);
  std::vector<arma::rowvec> e_group(c);
  std::vector<arma::mat> x_group(c);

  for (uword i=0; i<n; ++i)
    ++G_count[G[i]];

  for (uword i=0; i<c; ++i) {
    e_group[i].resize(G_count[i]);
    x_group[i].resize(G_count[i], k);
  }

  uword g;
  for (uword i=n; i>0; --i) {
    g = G[i-1];
    e_group[g][--G_count[g]] = e[i-1];
    x_group[g].row(G_count[g]) = X.row(i-1);
  }

  // combiner can also be simplified as omp_out += omp_in
  #pragma omp declare reduction( + : arma::mat : \
      std::transform(omp_in.begin(), omp_in.end(), \
                     omp_out.begin(), omp_out.begin(), \
                     std::plus<double>())) \
    initializer(omp_priv(omp_orig))


  arma::mat inner(k, k, arma::fill::zeros);
  arma::rowvec u;

  omp_set_num_threads(NUM_THREAD);
  #pragma omp parallel for private(u) reduction(+ : inner)
  for (uword i=0; i<c; ++i) {
    u = e_group[i] * x_group[i];
    inner += u.t() * u;
  }

  // need adjustment?
  // double dbl_n = static_cast<double>(n), dbl_k = static_cast<double>(k);
  // std::sqrt((dbl_n-dbl_k)/(dbl_n-1))

  arma::mat vcov = outer * inner * outer;
  arma::vec robust_stderr = arma::sqrt(vcov.diag());

  arma::vec t_statistics = beta / robust_stderr;

  arma::vec p_value = t_statistics;
  p_value.for_each([&df](double &x){ x=R::pt(std::abs(x), df, false, false); });

  return List::create(
    _["parameters"] = beta,
    _["vcov_matrix"] = vcov,
    _["robust_std.err"] = robust_stderr,
    _["p_value"] = p_value,
    _["t_statistics"] = t_statistics
  );
}