#include <RcppArmadillo.h>
#include <queue>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;
using arma::uword;

template <typename T>
arma::uvec top_i_heap(const arma::Col<T>& v, uword n) {
  typedef pair<T, uword> pt_;
  priority_queue<pt_, vector<pt_>, greater<pt_>> heap;

  pt_ p;
  for (uword i=0; i<v.size(); ++i) {
    if (heap.size()<n)
      heap.push(pt_(v[i], i));
    else {
      p = pt_(v[i], i);
      if (heap.top() < p) {
        heap.pop();
        heap.push(p);
      }
    }
  }

  arma::uvec res(n);
  for (uword i=0; i<n; ++i) {
    res[i] = heap.top().second + 1;
    heap.pop();
  }

  return res;
}


// [[Rcpp::export]]
arma::uvec top_index(SEXP x, uword n) {
  switch(TYPEOF(x)) {
    case INTSXP: return top_i_heap(as<arma::ivec>(x), n);
    case REALSXP: return top_i_heap(as<arma::vec>(x), n);
    case CPLXSXP: return top_i_heap(as<arma::cx_vec>(x), n);
    // case STRSXP: return top_i_heap(as<arma::Col<string>>(x), n);
    default: Rcpp::stop("unexpected data type");
  }
  return arma::uvec();
}
