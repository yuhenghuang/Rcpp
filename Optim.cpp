// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

#include <omp.h>
#include <RcppArmadillo.h>

#define OPTIM_USE_RCPP_ARMADILLO
#include <optim.hpp>


using namespace Rcpp;
using arma::uword;

double ackley_fn(const arma::vec &vals_inp, arma::vec *grad_out, void *opt_data) {
  const double x = vals_inp(0);
  const double y = vals_inp(1);
  const double pi = arma::datum::pi;

  double obj_val = -20 * std::exp(-0.2*std::sqrt(0.5*(x*x + y*y))) - std::exp(0.5*(std::cos(2*pi*x) + std::cos(2*pi*y))) + 22.718282L;

  return obj_val;
}

// [[Rcpp::export]]
int ackley_test() {
  arma::vec x = arma::ones(2, 1) + 1.0;

  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  bool success = optim::de(x, ackley_fn, nullptr);

  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  if (success) {
    Rcout << "de: Ackley test completed successfully.\n"
          << "elapsed time: " << elapsed_seconds.count() << "s\n";
  }
  else {
    Rcout << "de: Ackley test completed unsuccessfully." << std::endl;
  }

  Rcout << "\nde: solution to Ackley test:\n" << x << std::endl;

  return 0;
}

arma::mat inline sigm(const arma::mat &X) {
  return 1.0 / (1.0 + arma::exp(-X));
}

struct ll_data_t
{
  /* data */
  arma::vec Y;
  arma::mat X;
};

double ll_fn_whess(const arma::vec &vals_inp,
                   arma::vec *grad_out, 
                   arma::mat *hess_out,
                   void *opt_data) {
  // ...
  ll_data_t *objfn_data = reinterpret_cast<ll_data_t*>(opt_data);

  const arma::vec &Y = objfn_data->Y;
  const arma::mat &X = objfn_data->X;

  arma::vec mu = sigm(X * vals_inp);

  const double norm_term = static_cast<double>(Y.n_elem);

  const double obj_val = - arma::accu(Y % arma::log(mu) + (1.0-Y) % arma::log(1.0-mu)) / norm_term;

  if (grad_out) {
    *grad_out = X.t() * (mu-Y) / norm_term;
  }

  if (hess_out) {
    arma::mat S = arma::diagmat(mu % (1.0-mu));
    *hess_out = X.t() * S * X / norm_term;
  }

  return obj_val;
}

double ll_fn(const arma::vec &vals_inp, arma::vec *grad_out, void * opt_data) {
  return ll_fn_whess(vals_inp, grad_out, nullptr, opt_data);
}

// [[Rcpp::export]]
int logistic_test() {
  uword n_dim = 5, n_smpl = 4000;

  arma::mat X = arma::randn(n_smpl, n_dim);
  arma::vec theta_0 = 1.0 + 3.0 * arma::randu(n_dim, 1);

  arma::vec mu = sigm(X * theta_0);

  arma::vec Y(n_smpl);

  for (uword i=0; i<n_smpl; ++i)
    Y[i] = (R::unif_rand() < mu(i)) ? 1.0 : 0.0;
  
  ll_data_t opt_data;
  opt_data.Y = std::move(Y);
  opt_data.X = std::move(X);

  arma::vec x = arma::ones(n_dim, 1) + 1.0;

  optim::algo_settings_t settings;

  settings.gd_settings.method = 6;
  settings.gd_settings.par_step_size = 0.1;

  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  bool success = optim::gd(x,ll_fn,&opt_data,settings);

  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  if (success) {
    Rcout << "Adam: logit_reg test completed successfully.\n"
          << "elapsed time: " << elapsed_seconds.count() << "s\n";
  }
  else {
    Rcout << "Adam: logit_reg test completed unsuccessfully." << std::endl;
  }

  Rcout << "\nAdam: true values vs estimates:\n" << arma::join_rows(theta_0, x) << std::endl;

  x = arma::ones(n_dim, 1) + 1.0;

  start = std::chrono::system_clock::now();

  success = optim::newton(x,ll_fn_whess,&opt_data);

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;

  if (success) {
    Rcout << "newton: logit_reg test completed successfully.\n"
          << "elapsed time: " << elapsed_seconds.count() << "s\n";
  }
  else {
    Rcout << "newton: logit_reg test completed unsuccessfully." << std::endl;
  }

  Rcout << "\nnewton: true values vs estimates:\n" << arma::join_rows(theta_0, x) << std::endl;

  return 0;
}
