#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

class Tauchen {
  private:
    arma::vec converge(const int &n,
                       const arma::mat &P);
  public:
    Tauchen(const double rho,
            const double sigmu_u,
            const double b,
            const int m,
            const int n);
    arma::mat P;
    arma::vec state_values, stationary_distributions;
};

arma::vec Tauchen::converge(const int &n,
                            const arma::mat &P) {
  /*
  compute stationary distributions
  */
  double err = 1.;
  uword count = 0;
  arma::rowvec init = arma::ones(1, n) / n, recu;
  while (err>1e-7 && count<4000) {
    recu = init * P;

    if (count % 10 ==0)
      err = arma::abs(recu-init).max();

    init = recu;
    ++count;
  }

  return init.t();
}

Tauchen::Tauchen(const double rho,
                 const double sigmu_u,
                 const double b=0.0,
                 const int m=3,
                 const int n=7) {
  /*
  y_{t+1} = b + \rho y_t + u_{t+1}
  */

  double y_lim = m * sigmu_u / std::sqrt(1-rho*rho), d_in = y_lim*2/(n-1);
  Rcout << d_in << std::endl;
  arma::vec y = arma::linspace<arma::vec>(-y_lim, y_lim, n);

  arma::mat pi(n, n);
  for (uword i=0; i<n; ++i) {
    pi(i, 0) = arma::normcdf((y[0] - rho*y[i] + d_in/2) / sigmu_u);
    pi(i, n-1) = 1 - arma::normcdf((y[n-1] - rho*y[i] - d_in/2) / sigmu_u);
    for (uword j=1; j<n-1; ++j)
      pi(i, j) = arma::normcdf((y[j] - rho*y[i] + d_in/2) / sigmu_u) - arma::normcdf((y[j] - rho*y[i] - d_in/2) / sigmu_u);
  }

  pi.each_row( [](arma::rowvec& x){ x /= arma::sum(x); } );

  this->state_values = y;
  this->P = pi;
  this->stationary_distributions = converge(n, pi);
}

// [[Rcpp::export]]
List tauchen_test(const double rho,
                  const double sigmu_u,
                  const int m=3,
                  const int n=7) {
  // ...
  Tauchen tc(rho, sigmu_u, 0.0, m, n);
  return List::create(
    _["grid"] = tc.state_values,
    _["pi"] = tc.P,
    _["station"] = tc.stationary_distributions
  );
}