#include <mlpack/core.hpp>
#include <mlpack/core/data/split_data.hpp>

#include <mlpack/methods/ann/layer/layer.hpp>
#include <mlpack/methods/ann/ffn.hpp>

#include <ensmallen.hpp>
// ...
using namespace std;
using arma::uword;

arma::Row<size_t> get_labels(arma::mat &pred_out) {
  arma::Row<size_t> pred_lbls(pred_out.n_cols);
  for (uword i=0; i<pred_out.n_cols; ++i)
    pred_lbls[i] = pred_out.col(i).index_max() + 1;
  return pred_lbls;
}

int main() {
  // ...

  constexpr double RATIO = 1.0;

  constexpr int MAX_ITER = 0;

  constexpr double STEP_SIZE = 1.2e-3;

  constexpr int BATCH_SIZE = 50;

  arma::mat temp_dataset;

  mlpack::data::Load("../data/train.csv", temp_dataset, true);

  arma::mat dataset = temp_dataset.submat(0, 1, temp_dataset.n_rows-1, temp_dataset.n_cols-1);

  arma::mat train, valid;

  mlpack::data::Split()
}