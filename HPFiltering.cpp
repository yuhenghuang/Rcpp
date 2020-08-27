/*
question description
https://web.stanford.edu/~boyd/papers/pdf/l1_trend_filter.pdf

primal-dual interior-point method
http://www.stat.cmu.edu/~ryantibs/convexopt/lectures/primal-dual.pdf

algorithm of solving for banded-matrix
https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
*/




// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using arma::uword;


// Hyperparameters for search algorithm
static const double ALPHA = 0.7;
static const double BETA = 0.9;
static const double MU = 1.2;


// Parameters
struct Params {
  double lambda;
  arma::vec filter_d, filter_dd;
  arma::vec dd1, dd2, dd3, dy;

  Params(uword n, const arma::vec &y, double l): lambda(l) {
    dd1.resize(n);
    dd2.resize(n);
    dd3.resize(n);

    dd1.fill(6);
    dd2.fill(-4);
    dd3.fill(1);

    filter_dd = {1, -4, 6, -4, 1};

    filter_d = {1, -2, 1};
    dy.resize(n);
    dy.fill(0);
    for (uword i=0; i<n; ++i) 
      for (uword j=0; j<3; ++j)
        dy[i] += y[i+j] * filter_d[j];
  }
};


// solve the first row of equation (17)
arma::vec solve_banded_matrix(arma::vec a,
                              arma::vec b,
                              arma::vec &c,
                              arma::vec d,
                              arma::vec e,
                              arma::vec &f) {
  // ...
  uword n = f.n_elem;
  arma::vec res(n);
  double w;
  for (uword i=2; i<n; ++i) {
    w = a[i] / b[i-1];
    b[i] -= w * c[i-1];
    c[i] -= w * d[i-1];
    d[i] -= w * e[i-1];
    f[i] -= w * f[i-1];
  }

  for (uword i=1; i<n; ++i) {
    w = b[i] / c[i-1];
    c[i] -= w * d[i-1];
    d[i] -= w * e[i-1];
    f[i] -= w * f[i-1];
  }

  for (uword i=n-2; i>0; --i) {
    w = e[i-1] / d[i];
    d[i-1] -= w * c[i];
    f[i-1] -= w * f[i];
  }

  res[n-1] = f[n-1] / c[n-1];

  for (uword i=n-1; i>0; --i) {
    res[i-1] = (f[i-1] - d[i-1] * res[i]) / c[i-1];
  }

  return res;
}


// compute D D^T v
arma::vec DD_filter(const arma::vec &v) {
  uword n = v.n_elem;
  arma::vec res(n, arma::fill::zeros);

  arma::rowvec filter = {1, -4, 6, -4, 1};
  for (uword i=2; i<n-2; ++i) 
    for (uword j=0; j<5; ++j)
      res[i] += filter[j] * v[i-2+j];

  uword idx;
  for (uword i=0; i<2; ++i)
    for (uword j=0; j<5; ++j) {
      idx = i-2+j;
      if (idx>=0) {
        res[i] += filter[j] * v[idx];
        res[n-1-i] += filter[4-j] * v[n-1-idx];
      }
    }

  return res;
}


// for memory efficiency
void inline fn_f_inv(arma::vec &f, const double &lambda) {
  f.for_each( [&lambda](double x){ x = 1/(x-lambda); } );
}


// initialize s
double init_s(const arma::vec &mu1,
              const arma::vec &mu2,
              const arma::vec &dmu1,
              const arma::vec &dmu2) {
  // ...
  uword n = mu1.n_elem;
  double s = 1.0;
  for (uword i=0; i<n; ++i) {
    if (dmu1[i]<0)
      s = std::min(s, -mu1[i] / dmu1[i]);
    
    if (dmu2[i]<0)
      s = std::min(s, -mu2[i] / dmu2[i]);
  }

  return s * 0.999;
}


double compute_r_norm_plus(const arma::vec &v,
                           const arma::vec &mu1,
                           const arma::vec &mu2,
                           const arma::vec &dv,
                           const arma::vec &dmu1,
                           const arma::vec &dmu2,
                           double s,
                           double t,
                           Params *par) {
  // ...
  arma::vec v_plus = v + s * dv;
  arma::vec mu1_plus = mu1 + s * dmu1;
  arma::vec mu2_plus = mu2 + s * dmu2;

  arma::vec r_v = par->dy;
  // r_v -= DD_filter(v_plus);
  r_v -= arma::conv(v_plus, par->filter_dd, "same");
  r_v -= mu1_plus;
  r_v += mu2_plus;

  // modify in-place to avoid unnecessary copy
  mu1_plus %= (-par->lambda + v_plus);
  mu2_plus %= (-par->lambda - v_plus);

  mu1_plus += 1/t;
  mu2_plus += 1/t;

  return arma::dot(r_v.t(), r_v) + arma::dot(mu1_plus.t(), mu1_plus) + arma::dot(mu2_plus.t(), mu2_plus);
}


//' @export
// [[Rcpp::export]]
List hp_filter(const arma::vec &y,
               const double lambda) {
  // ...
  uword n = y.n_elem - 2;
  Params *par = new Params(n, y, lambda);

  arma::vec v(n, arma::fill::zeros);
  arma::vec mu1(n, arma::fill::ones);
  arma::vec mu2(n, arma::fill::ones);

  double eta = - arma::dot(mu1.t(), (-lambda+v)) - arma::dot(mu2.t(), (-lambda-v));


  int k = 1;
  double r_norm = 1.0, r_norm_plus = 1.0;

  double t, s;
  arma::vec J1_inv, J2_inv, f1_inv, f2_inv, v_lhs, v_rhs;
  arma::vec r_v, r_mu1, r_mu2;
  arma::vec dv, dmu1, dmu2;
  while (eta > 1e-6 || r_norm > 1e-6) {
    t = MU * n / std::pow(eta, k-1);

    // prepare parts for computing dv
    f1_inv = v;
    f2_inv = -v;
    fn_f_inv(f1_inv, lambda);
    fn_f_inv(f2_inv, lambda);

    J1_inv = f1_inv % mu1;
    J2_inv = f2_inv % mu2;

    v_rhs = par->dy;
    // v_rhs -= DD_filter(v);
    v_rhs -= arma::conv(v, par->filter_dd, "same");

    // save the temporary result for computing norm of r dual
    r_v = v_rhs;

    v_rhs += f1_inv / t;
    v_rhs -= f1_inv / t;

    v_lhs = par->dd1;
    v_lhs -= J1_inv;
    v_lhs += J2_inv;

    dv = solve_banded_matrix(par->dd3, par->dd2, v_lhs, par->dd2, par->dd3, v_rhs);

    dmu1 = - mu1 - f1_inv / t;
    dmu2 = - mu2 - f2_inv / t;
    dmu1 += J1_inv % dv;
    dmu2 += J2_inv % dv;


    // backtracking line search
    // search proper s
    // the selection of s ensures that updated mu1, mu2 are positive
    s = init_s(mu1, mu2, dmu1, dmu2);

    r_v -= mu1;
    r_v += mu2;
    r_mu1 = mu1 % (-lambda + v) + 1/t;
    r_mu2 = mu2 % (-lambda - v) + 1/t;

    r_norm = arma::dot(r_v.t(), r_v) + arma::dot(r_mu1.t(), r_mu1), arma::dot(r_mu2.t(), r_mu2);

    s /= BETA;
    do {
      s *= BETA;
      r_norm_plus = compute_r_norm_plus(v, mu1, mu2, dv, dmu1, dmu2, s, t, par);
    } while (r_norm_plus > (1-ALPHA*s) * r_norm);

    // update values
    v += s * dv;
    mu1 += s * dmu1;
    mu2 += s * dmu2;

    // compute error terms
    eta = - arma::dot(mu1.t(), (-lambda+v)) - arma::dot(mu2.t(), (-lambda-v));
    r_norm = std::sqrt(r_norm);
    ++k;
  }

  arma::vec x = y;
  x -= arma::conv(v, par->filter_d);

  return List::create(
    _["x"] = x
  );
}