#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
List em1(const arma::mat y,
         const arma::mat X,
         const int maxit = 10) {
  // inputs
  const int N = y.n_rows;
  const int K = X.n_cols;

  arma::mat beta(K, 1);
  beta.fill(0.);
  arma::mat eystar(N, 1);
  eystar.fill(0);

  for (int it=0; it<maxit; ++it) {
    // Reserved for next step
  }

  return List::create(
    _["N"] = N,
    _["K"] = K,
    _["beta"] = beta,
    _["eystar"] = eystar
  );
}

// [[Rcpp::export]]
List em2(const arma::mat y,
         const arma::mat X,
         const int maxit = 10) {
  // inputs
  const int N = y.n_rows;
  const int K = X.n_cols;

  arma::mat beta(K, 1);
  beta.fill(0.);
  arma::mat eystar(N, 1);
  eystar.fill(0);

  for (int it=0; it<maxit; ++it) {
    arma::mat mu = X * beta;
    // augmentation step
    for (int n=0; n<N; ++n)
      eystar(n, 0) = y(n, 0)==1 ? 1 : -1;

    // maximization step
    beta = (X.t() * X).i() * X.t() * eystar;
  }

  return List::create(
    _["N"] = N,
    _["K"] = K,
    _["beta"] = beta,
    _["eystar"] = eystar
  );
}

double inline f(double mu) {
  double val = R::dnorm4(-mu, 0, 1, false) /
               (1 - R::pnorm5(-mu, 0, 1, true, false));
  return mu + val;
}

double inline g(double mu) {
  double val = R::dnorm4(-mu, 0, 1, false) /
               R::pnorm5(-mu, 0, 1, true, false);
  return mu - val;
}

// [[Rcpp::export]]
List em3(const arma::mat y,
         const arma::mat X,
         const int maxit = 10) {
  // inputs
  const int N = y.n_rows;
  const int K = X.n_cols;

  arma::mat beta(K, 1);
  beta.fill(0.);
  arma::mat eystar(N, 1);
  eystar.fill(0);

  for (int it=0; it<maxit; ++it) {
    arma::mat mu = X * beta;
    // augmentation step
    for (int n=0; n<N; ++n)
      eystar(n, 0) = y(n, 0)==1 ? f(mu(n, 0)) : g(mu(n, 0));

    // maximization step
    beta = (X.t() * X).i() * X.t() * eystar;
  }

  return List::create(
    _["N"] = N,
    _["K"] = K,
    _["beta"] = beta,
    _["eystar"] = eystar
  );
}



// [[Rcpp::export]]
List em4(const arma::mat y,
         const arma::mat X,
         const int maxit = 10,
         const int nthr = 1) {
  // inputs
  const int N = y.n_rows;
  const int K = X.n_cols;
  omp_set_num_threads(nthr);

  arma::mat beta(K, 1);
  beta.fill(0.);
  arma::mat eystar(N, 1);
  eystar.fill(0);

  for (int it=0; it<maxit; ++it) {
    arma::mat mu = X * beta;
    // augmentation step
    #pragma omp parallel for
    for (int n=0; n<N; ++n)
      eystar(n, 0) = y(n, 0)==1 ? f(mu(n, 0)) : g(mu(n, 0));

    // maximization step
    beta = (X.t() * X).i() * X.t() * eystar;
  }

  return List::create(
    _["N"] = N,
    _["K"] = K,
    _["beta"] = beta,
    _["eystar"] = eystar
  );
}