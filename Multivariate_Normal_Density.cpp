#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

static const double log2pi = std::log(2. * M_PI);

// [[Rcpp::export]]
arma::vec dmvnorm_arma(const arma::mat &x,
                       const arma::rowvec &mu,
                       const arma::mat &sigma,
                       const bool logd = false) {
  
  /*
  x: (n, k)
  mu: (k,)
  sigma: (k, k)
  */
  // init
  // unsigned int32
  
  using arma::uword;

  const uword n = x.n_rows, k = x.n_cols;
  arma::vec out(n);
  const arma::mat rooti = arma::trimatu(arma::chol(sigma)).i();

  const double rootisum = arma::sum(log(rooti.diag())),
               constants = - double(k) / 2. * log2pi,
               other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i=0; i<n; ++i) {
    z = (x.row(i) - mu) * rooti;
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  return logd ? exp(out) : out;
}

void inplace_trimat_mul(arma::rowvec &x, const arma::mat &trimat) {
  const arma::uword n = trimat.n_cols;

  // be wary of unsigned int...
  for (arma::uword j=n; j-->0;) {
    double temp = 0.;
    for (arma::uword i=0; i<=j; ++i)
      temp += trimat(i, j) * x[i];
    x[j] = temp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma_fast(const arma::mat &x,
                            const arma::rowvec &mu,
                            const arma::mat &sigma,
                            const bool logd = false) {
  
  /*
  x: (n, k)
  mu: (k,)
  sigma: (k, k)
  */
  // init
  // unsigned int32
  
  using arma::uword;

  const uword n = x.n_rows, k = x.n_cols;
  arma::vec out(n);
  const arma::mat rooti = arma::trimatu(arma::chol(sigma)).i();

  const double rootisum = arma::sum(log(rooti.diag())),
               constants = - double(k) / 2. * log2pi,
               other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i=0; i<n; ++i) {
    z = x.row(i) - mu;
    inplace_trimat_mul(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  return logd ? exp(out) : out;
}


// [[Rcpp::export]]
arma::vec dmvnorm_arma_mc(const arma::mat &x,
                          const arma::rowvec &mu,
                          const arma::mat &sigma,
                          const bool logd = false,
                          const int cores = 1) {
  
  /*
  x: (n, k)
  mu: (k,)
  sigma: (k, k)
  */
  // init
  // unsigned int32
  omp_set_num_threads(cores);

  using arma::uword;

  const uword n = x.n_rows, k = x.n_cols;
  arma::vec out(n);
  const arma::mat rooti = arma::trimatu(arma::chol(sigma)).i();

  const double rootisum = arma::sum(log(rooti.diag())),
               constants = - double(k) / 2. * log2pi,
               other_terms = rootisum + constants;

  arma::rowvec z;
  #pragma omp parallel for schedule(static) private(z)
  for (uword i=0; i<n; ++i) {
    z = x.row(i) - mu;
    inplace_trimat_mul(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  return logd ? exp(out) : out;
}