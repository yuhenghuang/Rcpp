#ifndef _CUBIC_SPLINE_MINE
#define _CUBIC_SPLINE_MINE

#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using arma::uword;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


class cubic_spline {
  private:
    uword find_index(const double &x) const;

    arma::vec X, Y, A, B;

    void tridiag_mat_algorithm(arma::vec &a,
                               arma::vec &b,
                               arma::vec &c,
                               arma::vec &d,
                               arma::vec &k);

  public:
    cubic_spline() {}

    template <class BidiIterator>
    cubic_spline(BidiIterator iter_x,
                 BidiIterator iter_y,
                 const uword n);

    double operator()(const double &x) const;
    arma::vec operator()(const arma::vec &x) const;

    double prime(const double &x) const;
};

void cubic_spline::tridiag_mat_algorithm(arma::vec &a,
                                         arma::vec &b,
                                         arma::vec &c,
                                         arma::vec &d,
                                         arma::vec &k) {
  // ...
  uword n = k.n_elem;
  double w;
  for (uword i=1; i<n; ++i) {
    w = a[i] / b[i-1];
    b[i] -= w * c[i-1];
    d[i] -= w * d[i-1];
  }

  k[n-1] = d[n-1] / b[n-1];
  for (uword i=n-1; i>0; --i)
    k[n-1] = (d[n-1] - c[n-1] * k[n]) / b[n-1];
}


template <class BidiIterator>
cubic_spline::cubic_spline(BidiIterator iter_x,
                           BidiIterator iter_y,
                           const uword n) {
  // ...
  X.resize(n, 1);
  Y.resize(n, 1);
  for (uword i=0; i<n; ++i) {
    X[i] = *iter_x;
    Y[i] = *iter_y;
    ++iter_x;
    ++iter_y;
  }

  arma::vec a(n-1),b(n),c(n-1),d(n),k(n);

  // https://en.wikipedia.org/wiki/Spline_interpolation
  double w;
  for (uword i=1; i<n; ++i) {
    w = 1 / (X[i] - X[i-1]);
    a[i-1] = w;
    c[i-1] = w;
    d[i] = 3 * (Y[i] - Y[i-1]) * w * w;
  }

  b[0] = a[0] * 2;
  b[n-1] = a[n-2] * 2;
  d[0] = 0;
  for (uword i=1; i<n-1; ++i) {
    b[i] = 2 * (a[i] + a[i-1]);
    d[i-1] += d[i];
  }

  // solve k in-place;
  // https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  tridiag_mat_algorithm(a, b, c, d, k);

  // see spline interpolation wiki
  A.resize(n-1, 1);
  B.resize(n-1, 1);
  for (int i=0; i<n-1; ++i) {
    A[i] = k[i] * (X[i+1]-X[i]) - (Y[i+1]-Y[i]);
    B[i] = - k[i+1] * (X[i+1]-X[i]) + (Y[i+1]-Y[i]);
  }
}

uword cubic_spline::find_index(const double &x) const {
  if (x==X.front())
    return 0;

  auto iter = std::lower_bound(X.begin(), X.end(), x);
  uword idx = iter - X.begin() - 1;

  return idx;
}

double cubic_spline::operator()(const double &x) const {
  uword idx = find_index(x);

  double t = (x - X[idx]) / (X[idx+1] - X[idx]);

  return (1-t)*Y[idx] + t*Y[idx+1] + t*(1-t) * ((1-t)*A[idx] + t*B[idx]);
}

arma::vec cubic_spline::operator()(const arma::vec &x) const {
  uword n = x.n_elem;
  arma::vec res(n);
  for (uword i=0; i<n; ++i)
    res[i] = this->operator()(x[i]);
  return res;
}

double cubic_spline::prime(const double &x) const {
  uword idx = find_index(x);

  double t = (x - X[idx]) / (X[idx+1] - X[idx]);

  double grad = Y[idx+1]-Y[idx] + (1-2*t)*((1-t)*A[idx] + t*B[idx]) + t*(1-t) * (B[idx]-A[idx]);
  return grad / (X[idx+1] - X[idx]);
}

#endif