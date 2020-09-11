#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using arma::uword;

struct Kalman {
  arma::mat A, H, Q, R;
  double dt;

  Kalman(): dt(1.0) {
    // the state transition model
    A = arma::eye(6, 6);
    A(0, 2) = A(1, 3) = dt;
    A(2, 4) = A(3, 5) = dt;

    // the observation model
    H = arma::zeros(2, 6);
    H(0, 0) = H(1, 1) = dt;

    // the covariance of the process noise
    Q = arma::eye(6, 6);

    // the covariance of the observation noise
    R = arma::eye(2, 2) * 1000.0; 
  }

  arma::mat estimate(const arma::mat &Z) {
    // current state
    arma::vec x_est(6, arma::fill::zeros);
    // current covariance
    arma::mat p_est(6, 6, arma::fill::zeros);

    uword n = Z.n_rows, k = Z.n_cols;

    arma::mat Y(n, k);

    arma::mat p_prd, S, B, kg;
    arma::vec x_prd, z, y;

    for (uword i=0; i<n; ++i) {
      z = Z.row(i).t();

      // predict state and covariance
      x_prd = A * x_est;
      p_prd = A * p_est * A.t() + Q;

      // innovation (pre-fit residual) covariance
      B = H * p_prd.t();
      S = B * H.t() + R;
      kg = (arma::solve(S, B)).t();

      x_est = x_prd + kg * (z - H * x_prd);
      p_est = p_prd - kg * H * p_prd;

      y = H * x_est;
      Y.row(i) = y.t();
    }

    return Y;
  }
};


// [[Rcpp::export]]
arma::mat KalmanFilter(const arma::mat &Z) {
  Kalman K;
  return K.estimate(Z);
}
