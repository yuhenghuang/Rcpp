// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector convolveC(const NumericVector &a,
                        const NumericVector &b) {
  
  // init
  int N_a = a.size(), N_b = b.size(), N_ab = N_a + N_b - 1;

  NumericVector ab(N_ab);

  for (int i=0; i<N_a; ++i)
    for (int j=0; j<N_b; ++j) {
      ab[i+j] += a[i] * b[j];
    }

  return ab;
}

// [[Rcpp::export]]
NumericVector bootstrapC(const NumericVector &ds, const int B = 1000) {

  NumericMatrix boot_stat(B, 2);

  int n = ds.size();

  for (int i=0; i<B; ++i) {
    NumericVector smpl = ds[floor(runif(n, 0, n))];

    boot_stat(i, 0) = mean(smpl);
    boot_stat(i, 1) = sd(smpl);
  }
  return boot_stat;
}


// [[Rcpp::export]]
arma::mat rmvnorm(int n,
                  const arma::vec &mu,
                  const arma::mat &sigma) {
  
  // multivariate normal distribution
  unsigned int p = sigma.n_cols;

  NumericVector draw = rnorm(n*p);

  // reshape in place, copy_aux_mem = false
  arma::mat Z = arma::mat(draw.begin(), n, p, false, true);

  // repreat mu from (p,) to (n, p)
  // Cholesky decomposition of variance-covariance matrix
  // an algorithm to compute matrix version of sqrt()
  arma::mat Y = arma::repmat(mu, 1, n).t() + Z * arma::chol(sigma);

  return Y;
}


// [[Rcpp::export]]
List fastLM(const arma::mat &X,
            const arma::colvec &y) {
  
  // Linear model
  // # of obs. and regressors
  int n = X.n_rows, k = X.n_cols;

  arma::colvec coef = (X.t() * X).i() * X.t() * y;

  arma::colvec res = y - X*coef;

  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.) / (n-k);

  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(X.t()*X)));

  return List::create(
    _["coefficients"] = coef,
    _["stderr"] = std_err,
    _["df.residual"] = res
  );
}