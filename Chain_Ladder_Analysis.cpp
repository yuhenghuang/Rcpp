// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using arma::uword;

double getFactor(arma::mat &mTri,
                 uword idx) {
  // ...
  uword nRow = mTri.n_rows;
  const arma::mat &subMat = mTri.submat(0, idx, nRow-idx-2, idx+1);
  arma::rowvec colSum = arma::sum(subMat, 0);
  return colSum[1] / colSum[0];
}

arma::rowvec getFactors(arma::mat &mTri) {
  // ...
  uword nCol = mTri.n_cols;
  arma::rowvec factors(nCol-1);
  for (uword i=0; i<nCol-1; ++i)
    factors[i] = getFactor(mTri, i);

  return factors;
}

// [[Rcpp::export]]
arma::mat getChainSquare_cpp(arma::mat &mTri) {
  /*
  mTri needs fliping left and right to be an upper triangular matrix
  */
  uword nRow = mTri.n_rows, nCol = mTri.n_cols;
  arma::mat out(mTri.begin(), nRow, nCol, true);

  arma::rowvec factors = getFactors(mTri);
  // flip left and right
  arma::vec dAntiDiag = arma::diagvec(arma::fliplr(mTri));

  arma::rowvec prodVec;
  for (uword i=1; i<nRow; ++i) {
    // subvec is closed on both sides, as well as span...
    prodVec = arma::cumprod(factors.subvec(nCol-i-1, nCol-2));
    out(i, arma::span(nCol-i, nCol-1)) = dAntiDiag[i] * prodVec;
  }
  return out;
}