/*
  https://gallery.rcpp.org/articles/hierarchical-risk-parity/
*/

#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using arma::uword;

// [[Rcpp::export]]
arma::mat distance_mat_elementwise(const arma::mat &mat_corr) {
  // ...correlation matrix
  uword n = mat_corr.n_rows, j;
  arma::mat dist_mat(n, n);

  // omp_set_num_threads(2);

  double d;
  #pragma omp parallel for private(j, d) collapse(2)
  for (uword i=1; i<n; ++i)
    for (j=0; j<i; ++j) {
      d = std::pow(0.5*(1-mat_corr(i, j)), 0.5);
      dist_mat(i, j) = d;
      dist_mat(j, i) = d;
    }

  return dist_mat;
}

// [[Rcpp::export]]
arma::mat distance_mat_rowwise(const arma::mat &mat_corr) {
  // ...
  uword n = mat_corr.n_rows, i, j, k;
  double temp_sum;
  arma::mat dist_mat(n, n);

  #pragma omp parallel for private(temp_sum, j, k)
  for (i=1; i<n; ++i)
    for (j=0; j<i; ++j) {
      temp_sum = 0;
      for (k=0; k<n; ++k)
        temp_sum += std::pow(mat_corr(k, i) - mat_corr(k, j), 2);
      temp_sum = std::pow(temp_sum, 0.5);
      dist_mat(i, j) = temp_sum;
      dist_mat(j, i) = temp_sum;
    }

  return dist_mat;
}

// [[Rcpp::export]]
arma::mat cluster_matrix(arma::mat mat_corr) {
  // ...
  uword dim = mat_corr.n_rows;
  double max_element = mat_corr.max(), min_dist;

  uword idx_row=0, idx_col=0;
  arma::mat temp_cluster_mat;
  arma::rowvec temp_cluster_row;
  arma::colvec temp_cluster_col;

  arma::mat cluster_mat(dim-1, 4);
  arma::colvec cluster_idx((dim-1)*2);

  #pragma omp parallel for
  for (uword i=0; i<dim; ++i)
    mat_corr(i, i) = max_element*2;

  #pragma omp parallel for
  for (uword i=0; i<(dim-1)*2; ++i)
    cluster_idx[i] = i+1;

  for (uword i=0; i<dim-1; ++i) {
    min_dist = mat_corr.min(idx_row, idx_col);

    cluster_mat(i, 0) = cluster_idx[idx_row];
    cluster_mat(i, 1) = cluster_idx[idx_col];
    cluster_mat(i, 2) = min_dist;
    cluster_mat(i, 3) = (cluster_mat(i,0)<=dim ? 1 : 0) + 
                           (cluster_mat(i,1)<=dim ? 1 : 0);

    cluster_idx.shed_row(idx_row);
    cluster_idx.shed_col(idx_col);

    temp_cluster_mat = arma::join_rows(mat_corr.col(idx_row), mat_corr.col(idx_col));

    temp_cluster_col = min(temp_cluster_mat, 1);
    temp_cluster_row = temp_cluster_col.t();
    // add one elements at the end of rowvec, 
    // set it to be max_element, as it will be diagonal element.
    temp_cluster_row.insert_cols(temp_cluster_col.n_elem, 1);
    temp_cluster_row[temp_cluster_row.n_elem-1] = max_element;

    mat_corr = arma::join_rows(mat_corr, temp_cluster_col);
    mat_corr = arma::join_cols(mat_corr, temp_cluster_row);

    mat_corr.shed_row(idx_row);
    mat_corr.shed_row(idx_col);
    mat_corr.shed_col(idx_row);
    mat_corr.shed_col(idx_col);
  }
  return cluster_mat;
}

arma::mat flat_cluster(uword index, 
                       uword &threshold, 
                       arma::mat &cluster_mat) {
  // ...
  arma::mat temp_idx, temp_idx_a, temp_idx_b;

  if (index<=threshold) {
    temp_idx.set_size(1, 1);
    temp_idx(0, 0) = index;
    return temp_idx;
  }

  temp_idx_a = flat_cluster(cluster_mat(index-threshold-1, 1), threshold, cluster_mat);
  temp_idx_b = flat_cluster(cluster_mat(index-threshold-1, 0), threshold, cluster_mat);
  temp_idx = arma::join_rows(temp_idx_a, temp_idx_b);
  return temp_idx;
}


// [[Rcpp::export]]
arma::mat cluster_index(arma::mat &cluster_mat, 
                        arma::mat &mat_cov) {
  // ...
  uword num_asset = mat_cov.n_rows, nrow_cluster = cluster_mat.n_rows;

  arma::mat asset_index = flat_cluster(nrow_cluster+num_asset, num_asset, cluster_mat);
  return asset_index;
}

// [[Rcpp::export]]
arma::mat quasi_diag(arma::mat &mat_cov, 
                     arma::mat &asset_index) {
  // ...
  uword num_asset = mat_cov.n_rows, index_asset;

  arma::mat inter_mat(num_asset, num_asset);
  arma::mat quasi_diagmat(num_asset, num_asset);

  #pragma omp parallel for private(index_asset)
  for (uword i=0; i<num_asset; ++i) {
    index_asset = asset_index(0, i) - 1;
    inter_mat.col(i) = mat_cov.col(index_asset);
  }

  #pragma omp parallel for private(index_asset)
  for (uword i=0; i<num_asset; ++i) {
    index_asset = asset_index(0, i) - 1;
    quasi_diagmat.row(i) = inter_mat.row(index_asset);
  }

  return quasi_diagmat;
}


void bisect_weight_allocation(arma::mat &mat_weight,
                              arma::mat &mat_cov,
                              uword idx_start,
                              uword idx_end) {
  // ...
  arma::colvec wi_upper, wi_lower;

  arma::mat mat_cov_upper, mat_cov_lower;

  uword idx_mid;

  double scale_upper, scale_lower;

  if (idx_start!=idx_end) {
    idx_mid = (idx_start+idx_end) / 2;

    mat_cov_upper = mat_cov.submat(idx_start, idx_start, idx_mid, idx_mid);
    mat_cov_lower = mat_cov.submat(idx_mid+1, idx_mid+1, idx_end, idx_end);

    wi_upper = mat_cov_upper.diag();
    wi_lower = mat_cov_upper.diag();

    scale_upper = as_scalar(wi_upper.t() * mat_cov_upper * wi_upper);
    scale_lower = as_scalar(wi_lower.t() * mat_cov_lower * wi_lower);

    mat_weight.submat(0, idx_start, 0, idx_mid) = mat_weight.submat(0, idx_start, 0, idx_mid) * 
                                                  (scale_lower / (scale_upper + scale_lower));
    mat_weight.submat(0, idx_mid+1, 0, idx_end) = mat_weight.submat(0, idx_mid+1, 0, idx_end) *
                                                  (scale_upper / (scale_upper + scale_lower));

    #pragma omp task shared(mat_weight, mat_cov) firstprivate(idx_start, idx_mid) 
    {
      bisect_weight_allocation(mat_weight, mat_cov, idx_start, idx_mid);
    }

    #pragma omp task shared(mat_weight, mat_cov) firstprivate(idx_mid, idx_end) 
    {
      bisect_weight_allocation(mat_weight, mat_cov, idx_mid+1, idx_end);
    }
  }
}

// [[Rcpp::export]]
arma::mat weight_allocation(arma::mat &quasi_diagmat,
                            arma::mat &asset_index) {
  // ...
  uword num_asset = quasi_diagmat.n_rows;

  arma::mat mat_cov(quasi_diagmat.begin(), num_asset, num_asset, false);
  arma::mat mat_weight_temp(1, num_asset, arma::fill::ones);
  arma::mat mat_weight(1, num_asset, arma::fill::ones);

  omp_set_nested(0);

  #pragma omp parallel
  {
    #pragma omp single
    {
      bisect_weight_allocation(mat_weight_temp, mat_cov, 0, num_asset-1);
    }
  }

  for (uword i=0; i<num_asset; ++i)
    mat_weight[0, i] = mat_weight_temp[0, asset_index[0, i]-1];

  return mat_weight;
}