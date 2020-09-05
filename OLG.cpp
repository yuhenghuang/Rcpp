#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

struct Params {
  double alpha, delta, beta;
  uword j_max;

  double k, w, r;
  arma::vec c, a;

  Params(double al, double d, double b, uword j): alpha(al), delta(d), beta(b), j_max(j) {
    c.resize(j_max);
    a.resize(j_max+1);
    a[0] = 0.0;
  }

  void update(double rate) {
    r = rate;
    k = std::pow((r+delta)/alpha, 1.0/(alpha-1.0)) * j_max;
    w = std::pow((1.0-alpha)*k/j_max, alpha);
  }
};

void OLG_inner(Params *par) {

  double discount = 0.0, c_agg = 0.0;
  arma::vec &c = par->c;
  arma::vec &a = par->a;

  double multiplier_c = 1.0, multiplier_i = 1.0;
  for (uword j=0; j<par->j_max; ++j) {
    c_agg += par->w / multiplier_i;
    discount += multiplier_c / multiplier_i;
    c[j] = multiplier_c;

    multiplier_c *= par->beta * (1+par->r);
    multiplier_i *= 1+par->r;
  }

  c *= c_agg / discount;
  for (uword j=0; j<par->j_max; ++j)
    a[j+1] = (1+par->r) * a[j] + par->w - c[j];
}


// [[Rcpp::export]]
List OLG_model(double alpha,
               double delta,
               double beta,
               uword j_max) {
  // ...
  Params *par = new Params(alpha, delta, beta, j_max);

  uword count = 0;
  double err = 1.0;
  double r_up = 1.0, r_low = -0.1;

  double kprime;
  while (err>1e-7 && count<100) {
    par->update((r_up+r_low)*0.5);

    OLG_inner(par);
    kprime = arma::sum(par->a);

    if (kprime < par->k)
      r_low = par->r;
    else
      r_up = par->r;

    err = std::abs(kprime - par->k);
    ++count;
  }

  return List::create(
    _["k"] = par->k,
    _["a"] = par->a,
    _["r"] = par->r,
    _["c"] = par->c
  );
}
