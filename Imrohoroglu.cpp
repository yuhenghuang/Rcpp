#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

arma::cube value_iter(arma::mat &v,
                      const arma::cube &c,
                      const arma::mat &pi,
                      const double &beta,
                      const double &delta) {
  // ...
  Rcout << "start value function iteration..." << std::endl;
  omp_set_num_threads(2);

  uword n=0, num_grids=v.n_rows, i;
  double err = 1.;
  arma::mat vn, expc;
  arma::cube u_base = arma::pow(c, 1-delta) / (1-delta), u;
  while (err>1e-7 && n<4000) {
    vn = v;
    // (4, 4) * (4, 301)
    expc = pi * v.t();
    u = u_base;
    // index inclusive at end in submat
    #pragma omp parallel for
    for (i=0; i<num_grids; ++i)
      u.slice(i) += beta * arma::repmat(expc.col(i).t(), num_grids, 1);

    v = arma::max(u, 2);

    if (n % 10 == 0)
      err = arma::abs(v-vn).max();

    if (n % 200 == 0)
      Rcout << "value function iteration at step " << n << " error is "<< err << std::endl;

    ++n;
  }

  Rcout << "value function iteration converged at " << n << " times" << std::endl;

  return u;
}

void demo_iter(arma::mat &demo,
               const arma::umat &plc,
               const arma::mat &pi) {
  // ...
  Rcout << "start demographic loops..." << std::endl;
  omp_set_num_threads(2);

  uword n=0, num_grids=demo.n_rows, choice;
  double err = 1.;
  arma::mat demon;
  while (err>1e-5 && n<4000) {
    demon = arma::zeros(num_grids, 4);

    #pragma omp parallel for private(choice)
    for (uword i=0; i<num_grids; ++i)
      for (uword j=0; j<4; ++j) {
        choice = plc(i, j);
        demon.row(choice) += demo(i, j) * pi.row(j);
      }

    if (n % 10 ==0)
      err = arma::abs(demon - demo).max();

    if (n % 200 == 0)
      Rcout << "demographic loops at step " << n << " error is "<< err << std::endl;

    demo = demon;
    ++n;
  }

  Rcout << "demographic loops converged at " << n << " times" << std::endl;
}


// [[Rcpp::export]]
List imrohoroglu(const double &beta,
                 const double &delta,
                 const double &theta,
                 const uword num_grids = 301) {
  // ...
  arma::mat bc_tran = {{0.9375, 0.0625}, {0.0625, 0.9375}};
  bc_tran = arma::repmat(bc_tran, 2, 2);
  bc_tran.swap_rows(1, 2);
  bc_tran.swap_cols(1, 2);

  arma::mat pi_good = {{0.975, 0.025}, {0.6, 0.4}};
  arma::mat pi_bad = {{145.0/154.0, 9.0/154.0}, {3.0/7.0, 4.0/7.0}};

  arma::mat pi = bc_tran % arma::repmat(arma::join_rows(pi_good, pi_bad), 2, 1);

  // transition matrix verified
  // Rcout << pi << std::endl;

  // value matrix
  // asset grid
  // income grid
  arma::mat v(num_grids, 4, arma::fill::ones);
  arma::vec a = arma::linspace<arma::vec>(0.0, 8.0, num_grids);
  arma::mat y = arma::repmat(arma::join_rows(arma::ones(num_grids, 1), arma::ones(num_grids, 1) * theta), 1, 2);

  // consumption
  arma::cube c(num_grids, 4, num_grids);
  arma::mat income = y + arma::repmat(a, 1, 4);
  c.each_slice() = income;

  // index inclusive at end in submat
  for (uword i=0; i<num_grids; ++i)
    c.slice(i).for_each( [&a, &i](double &x){ x -= a[i];} );

  // clip consumptions
  double c_min = 1e-8, c_max = 10;
  c.for_each( [&c_min, &c_max](double &x) { if (x<c_min) x=c_min; else if(x>c_max) x=c_max;});

  arma::cube u = value_iter(v, c, pi, beta, delta);

  arma::umat plc = arma::index_max(u, 2);

  arma::mat demo(num_grids, 4, arma::fill::ones);
  demo /= static_cast<double>(num_grids) * 4.0;

  demo_iter(demo, plc, pi);

  return List::create(
    _["a_grid"] = a,
    _["value"] = v,
    _["policy"] = plc,
    _["demo"] = demo
  );
}