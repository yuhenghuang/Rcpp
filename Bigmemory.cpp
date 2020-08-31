#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

#include <bigmemory/BigMatrix.h>

using namespace Rcpp;
using arma::uword;

template <typename T>
arma::Row<T> big_col_sum(const arma::Mat<T> &big_mat,
                         const arma::uvec &subset_cols) {
  // ...
  return arma::sum(big_mat.cols(subset_cols), 0);
}


// [[Rcpp::export]]
arma::rowvec big_col_sum(SEXP ptr_big_mat,
                         arma::uvec subset_cols) {
  // ...

  XPtr<BigMatrix> xp_mat(ptr_big_mat);

  if (any(subset_cols < 1) || any(subset_cols > xp_mat->ncol())) {
    throw std::out_of_range("Some of requested columns are outside of the matrix!");
  }

  subset_cols -= 1;

  unsigned int type = xp_mat->matrix_type();

  if (type==4) {
    arma::Row<int> col_sum = big_col_sum(
      arma::Mat<int>((int *)xp_mat->matrix(), xp_mat->nrow(), xp_mat->ncol(), false),
      subset_cols
    );
    return arma::conv_to<arma::rowvec>::from(col_sum);
  }
  else if (type==8) {
    arma::rowvec col_sum = big_col_sum(
      arma::mat((double *)xp_mat->matrix(), xp_mat->nrow(), xp_mat->ncol(), false),
      subset_cols
    );
    return col_sum;
  }
  else {
    throw Rcpp::exception("Undefined type encountered...");
  }
}
