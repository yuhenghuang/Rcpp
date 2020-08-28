#include <algorithm>

#include <RcppArmadillo.h>
#include <omp.h>
#include <Tauchen.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using arma::uword;

struct Params {
  const uword num_a_grid = 601,
              num_l_grid = 7;

  const double beta = 0.96,
               alpha = 0.36,
               delta = 0.08;
          
  double rho, mu, sigma;

  const double a_max = 40.0,
               b = 3.0;

  double l_min, phi, r, w, L, k;

  arma::vec a, l;

  arma::mat pi;

  Params(double r, double m, double s): rho(r), mu(m), sigma(s) {
    Tauchen x(rho, sigma*std::sqrt(1-rho*rho), 0, 3, num_l_grid);
    l = arma::exp(x.state_values);
    L = arma::sum(l % x.stationary_distributions);

    l_min = l[0];
    pi = x.P;
  }

  // update variables given rate
  void update(double rate) {
    r = rate;
    k = std::pow((r+delta) / (alpha * std::pow(L, 1.0-alpha)), 1.0/(alpha-1.0));
    w = (1.0-alpha) * std::pow(k, alpha) * std::pow(L, -alpha);

    phi = r>0 ? std::min(b, w*l_min/r) : b;

    a = arma::pow(arma::linspace(0, 1, num_a_grid), 2.0) * (phi+a_max) - phi;
  }
};

void plc_inner(Params *par,
               const arma::mat &c,
               const arma::mat &aprime,
               arma::mat &plc_c) {
  // ...
  uword sup;
  arma::vec temp;

  #pragma omp parallel for private(sup, temp)
  for (uword j=0; j<par->num_l_grid; ++j) {
    sup = std::lower_bound(par->a.begin(), par->a.end(), aprime(0, j)) - par->a.begin();

    arma::interp1(aprime.col(j), c.col(j), par->a.rows(sup, par->num_a_grid-1), temp, "linear", c(par->num_a_grid-1, j));
    plc_c.col(j).rows(sup, par->num_a_grid-1) = temp;

    if (sup>0) {
      plc_c.col(j).rows(0, sup-1) = (1.0+par->r) * par->a.rows(0, sup-1);
      plc_c.col(j).rows(0, sup-1) += par->phi + par->w * par->l[j];
    }
  }
}

arma::mat plc_iter(Params *par) {
  arma::mat plc_c(par->num_a_grid, par->num_l_grid, arma::fill::ones);
  plc_c *= 0.5 * par->a_max;

  uword count = 0;
  double err = 1.0;

  /*
  In this algorithm
  aprime is actually a_t
  par->a is a_t+1
  */
  arma::mat c, aprime, plc_c_orig;

  while (err > 1e-5 && count < 1000) {
    plc_c_orig = plc_c;
    c = arma::pow(plc_c, -par->mu);
    c *= par->pi.t();
    c *= par->beta * (1.0+par->r);
    c.for_each( [par](double &x){ x = std::pow(x, -1.0/par->mu); } );

    aprime = c;
    aprime.each_col() += par->a;
    aprime.each_row() -= par->w * par->l.t();
    aprime /= 1.0+par->r;

    plc_inner(par, c, aprime, plc_c);

    if (count % 5 ==0)
      err = arma::abs(plc_c - plc_c_orig).max();

    ++count;
  }

  // from here actually it is policy for aprime
  plc_c.for_each( [](double &x){ x=-x;} );
  plc_c.each_col() += (1.0+par->r) * par->a;
  plc_c.each_row() += par->w * par->l.t();
  return plc_c;
}

arma::mat demo_iter(Params *par,
                    const arma::mat &plc) {
  // ...
  arma::mat demo(par->num_a_grid, par->num_l_grid, arma::fill::ones);
  demo /= (par->num_a_grid * par->num_l_grid);

  arma::vec itv = par->a.rows(1, par->num_a_grid-1) - par->a.rows(0, par->num_a_grid-2);

  uword count = 0, idx;
  double err = 1.0;
  arma::rowvec temp;
  arma::mat demo_orig;
  while (err > 1e-6 && count < 1000) {
    demo_orig = demo;
    demo.fill(arma::fill::zeros);

    for (uword j=0; j<par->num_l_grid; ++j)

      #pragma omp parallel for private(idx, temp)
      for (uword i=0; i<par->num_a_grid; ++i) {
        temp = demo_orig(i, j) * par->pi.row(j);
        if (plc(i, j)==par->a[0])
          demo.row(0) += temp;
        else {
          idx = std::lower_bound(par->a.begin(), par->a.end(), plc(i, j)) - par->a.begin();
          demo.row(idx-1) += temp * ((par->a[idx] - plc(i,j)) / itv[idx-1]);
          demo.row(idx) += temp * ((plc(i,j) - par->a[idx-1]) / itv[idx-1]);
        }
      }
    
    if (count % 10 == 0)
      err = arma::abs(demo - demo_orig).max();

    ++count;
  }

  return demo;
}

//' @export
// [[Rcpp::export]]
List Aiyagari(double rho,
              double mu,
              double sigma) {
  // ...

  // set num of cores
  omp_set_num_threads(2);

  Params *par = new Params(rho, mu, sigma);

  double err = 1.0;
  uword count = 0;

  double r_low = -0.04, r_up = (1-par->beta) / par->beta, knew;
  arma::mat plc, demo;
  while (err>1e-4 && count<20) {
    par->update((r_low + r_up) / 2.0);

    plc = plc_iter(par);
    plc.for_each( [par](double &x){ if (x>par->a_max) x=par->a_max; else if (x<-par->phi) x=-par->phi; } );

    demo = demo_iter(par, plc);

    knew = arma::accu(demo % plc);

    err = std::abs(knew - par->k);

    if (knew < par->k)
      r_low = par->r;
    else 
      r_up = par->r;

    ++count;
  }

  return List::create(
    _["k"] = par->k,
    _["knew"] = knew,
    _["r"] = par->r,
    _["w"] = par->w,
    _["demo"] = demo,
    _["policy"] = plc
  );
}