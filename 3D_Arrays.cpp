#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

// [[Rcpp::export]]
arma::cube cube_omp(arma::cube &a,
                    const int cores = 1) {
  // ...
  uword xdim = a.n_rows, ydim = a.n_cols, tdim = a.n_slices;
  arma::cube res = a;

  omp_set_num_threads(cores);

  #pragma omp parallel for
  for (uword t=0; t<tdim; ++t) {
    arma::mat &temp_mat = a.slice(t);

    for (uword x=2; x<xdim-2; ++x) {
      const arma::mat &temp_row_sub = temp_mat.rows(x-2, x+2);

      for (uword y=2; y<ydim-2; ++y)
        res(x, y, t) = arma::accu(temp_row_sub.cols(y-2, y+2));
    }
  }
  return res;
}


// [[Rcpp::export]]
arma::cube cube_conv(arma::cube &a,
                     const int cores = 1) {
  // ...
  uword xdim = a.n_rows, ydim = a.n_cols, tdim = a.n_slices;
  arma::cube res(xdim, ydim, tdim);
  arma::mat filter(5, 5);
  filter.fill(1);

  omp_set_num_threads(cores);
  #pragma omp parallel for
  for (uword t=0; t<tdim; ++t) {
    res.slice(t) = arma::conv2(a.slice(t), filter, "same");
  }

  for (uword i=0; i<2; ++i) {
    res.row(i) = a.row(i);
    res.col(i) = a.col(i);
    res.row(xdim-i-1) = a.row(xdim-i-1);
    res.col(ydim-i-1) = a.col(ydim-i-1);
  }
  return res;
}