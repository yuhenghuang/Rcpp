#include <RcppArmadillo.h>
#include <vector>
#include <optim.hpp>
#include <omp.h>
#include <CubicSpline.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

template<typename T>
void clip(T &inp, const double &a_min, const double &a_max) {
  inp.for_each([&a_min, &a_max](double &x){if (x<a_min) x=a_min; else if(x>a_max) x=a_max;});
}

struct Params {
  // fixed parameters
  const uword num_grid_k_uneven = 130,
              num_grid_kbar_uneven = 6,
              num_grid_k_even = 600,
              num_grid_kbar_even = 100;

  const double k_min = 0.0, k_max = 26.0,
               kbar_min = 0.0, kbar_max = 13.0;

  const double z_good = 1.01,
               z_bad = 0.99;

  // concentration of capital grid points near 0
  const double skewness = 1.5;

  double beta, delta, sigma, alpha;

  // k, kbar grid points
  arma::vec k_uneven, k_even, kbar_uneven, kbar_even;

  // trasition matrix, hyperparameters
  arma::mat pi, H;

  // income
  arma::cube y_value, y_policy;

  // kbar prime grid points (2, num_grids...)
  arma::mat kbar_prime_value, kbar_prime_policy;

  Params(double &b, 
         double &d,
         double &s,
         double &a): beta(b), delta(d), sigma(s), alpha(a) {
    // grid_points
    k_uneven = arma::pow(arma::linspace<arma::vec>(0, 1, num_grid_k_uneven), skewness) * (k_max-k_min) + k_min;
    k_even = arma::linspace<arma::vec>(k_min, k_max, num_grid_k_even);

    kbar_uneven = arma::linspace<arma::vec>(kbar_min, kbar_max, num_grid_kbar_uneven);
    kbar_even = arma::linspace<arma::vec>(kbar_min, kbar_max, num_grid_kbar_even);

    // business cycle
    arma::mat bc = {{0.875, 0.125}, {0.125, 0.875}};

    arma::mat pi_gg = {{35.0/36.0, 1.0/36.0}, {2.0/3.0, 1.0/3.0}};
    arma::mat pi_bb = {{14.0/15.0, 2.0/30.0}, {0.6, 0.4}};

    arma::mat pi_gb = {{11.0/12.0, 1.0/12.0}, {0.5, 0.5}};
    arma::mat pi_bg = {{59.0/60.0, 1.0/60.0}, {0.75, 0.25}};

    arma::mat bc_tile = arma::repmat(bc, 2, 2);
    bc_tile.swap_rows(1, 2);
    bc_tile.swap_cols(1, 2);

    // ge, gu, be, bu
    pi = bc_tile % arma::join_cols(arma::join_rows(pi_gg, pi_gb), arma::join_rows(pi_bg, pi_bb));

    // hyper parameters
    H.resize(2, 2);
    setH(0.0, 0.95, 0.0, 0.94);

    // income
    y_value.resize(num_grid_k_uneven, num_grid_kbar_uneven, 4);
    y_policy.resize(num_grid_k_even, num_grid_kbar_even, 4);

    // kbar_prime;
    kbar_prime_value.resize(2, num_grid_kbar_uneven);
    kbar_prime_policy.resize(2, num_grid_kbar_even);

    // income
    double z, epsilon;
    arma::mat capital, wage;
    for (uword k=0; k<4; ++k) {
      z = (k<2) ? z_good : z_bad;
      epsilon = (k%2==0) ? 1 : 0;
      
      capital = arma::repmat(alpha*z*arma::pow(kbar_uneven.t(), alpha-1) +1-delta, num_grid_k_uneven, 1) % arma::repmat(k_uneven, 1, num_grid_kbar_uneven);
      wage = arma::repmat((1-alpha)*z*arma::pow(kbar_uneven.t() * epsilon, alpha), num_grid_k_uneven, 1);

      y_value.slice(k) = capital + wage;

      capital = arma::repmat(alpha*z*arma::pow(kbar_even.t(), alpha-1) +1-delta, num_grid_k_even, 1) % arma::repmat(k_even, 1, num_grid_kbar_even);
      wage = arma::repmat((1-alpha)*z*arma::pow(kbar_even.t() * epsilon, alpha), num_grid_k_even, 1);

      y_policy.slice(k) = capital + wage;
    }
  }

  void setH(double a0, double a1, double b0, double b1) {
    H(0, 0) = a0;
    H(0, 1) = a1;
    H(1, 0) = b0;
    H(1, 1) = b1;

    kbar_prime_value.row(0) = arma::exp(a0 + arma::log(kbar_uneven.t()) * a1);
    kbar_prime_value.row(1) = arma::exp(b0 + arma::log(kbar_uneven.t()) * b1);
    kbar_prime_policy.row(0) = arma::exp(a0 + arma::log(kbar_even.t()) * a1);
    kbar_prime_policy.row(1) = arma::exp(b0 + arma::log(kbar_even.t()) * b1);

    clip(kbar_prime_value, kbar_min, kbar_max);
    clip(kbar_prime_policy, kbar_min, kbar_max);
  }
};

struct StateVar {
  double beta, sigma, y, k_min, k_max;
  arma::rowvec trans;
  std::vector<cubic_spline> css;

  StateVar(double &b,
           double &s,
           const double &l,
           const double &r): beta(b), sigma(s), k_min(l), k_max(r) {
    // ...
    css.resize(4);
  }
};


double inline utility(const double &c,
                      const double &sigma) {
  // ...
  if (sigma==1.0)
    return std::log(c);
  else
    return (std::pow(c, 1-sigma) - 1) / (1-sigma);
}

double inline u_prime(const double &c,
                      const double &sigma) {
  // first derivative of utility function
  return std::pow(c, -sigma);
}

arma::mat inline utility(const arma::mat &c,
                         const double &sigma) {
  // ...
  if (sigma==1.0)
    return arma::log(c);
  else
    return (arma::pow(c, 1-sigma) - 1) / (1-sigma);
}

double val_fn(const arma::vec &vals_inp,
              arma::vec *grad_out,
              void *opt_data) {
  // minimization problem
  StateVar *s = reinterpret_cast<StateVar*>(opt_data);

  double kprime = vals_inp[0];
  double c = s->y - kprime;

  double obj_val = -utility(c, s->sigma),
         grad = -u_prime(c, s->sigma);

  double coef;
  for (uword l=0; l<4; ++l) {
    coef = s->beta * s->trans[l];
    obj_val -= coef * s->css[l](kprime);
    grad -= coef * s->css[l].prime(kprime);
  }

  if (grad_out) {
    (*grad_out)[0] = grad;
  }

  return obj_val;
}

arma::vec constr_fn(const arma::vec &vals_inp,
                    arma::mat *jacob_out,
                    void *opt_data) {
  // ...
  StateVar *s = reinterpret_cast<StateVar*>(opt_data);
  double kprime = vals_inp[0];
  arma::vec constr_vals(3);
  constr_vals[0] = s->k_min - kprime;
  constr_vals[1] = kprime - s->k_max;
  constr_vals[2] = kprime - s->y + 1e-6;

  if (jacob_out) {
    jacob_out->set_size(3, 1);

    (*jacob_out)(0, 0) = -1.0;
    (*jacob_out)(1, 0) = 1.0;
    (*jacob_out)(2, 0) = 1.0;
  }

  return constr_vals;
}

double bisection(StateVar *s) {
  double l = s->k_min, r = std::min(s->k_max, s->y) - 1e-6;
  double c, m, grad_u, grad_e;
  for (uword i=0; i<20; ++i) {
    m = (l+r) / 2;

    c = s->y - m;
    grad_u = u_prime(c, s->sigma);

    grad_e = 0;
    for (uword l=0; l<4; ++l)
      grad_e += s->beta * s->trans[l] * s->css[l].prime(m);

    if (std::abs(grad_e-grad_u)<1e-5)
      break;

    if (grad_e < grad_u)
      l = m;
    else
      r = m;
  }

  return m;
}

void value_to_policy(Params *par,
                     const arma::cube &value,
                     arma::cube &policy) {
  // ...
  StateVar *state = new StateVar(par->beta, par->sigma, par->k_min, par->k_max);

  bool success;
  arma::vec inp = {0.};

  // the shape is not wrong
  arma::cube value_prime(par->num_grid_k_uneven, par->num_grid_kbar_even, 4);

  for (uword k=0; k<4; ++k) {
    const arma::rowvec &kbar_prime = par->kbar_prime_policy.row(k/2);

    // only update twice (good and bad times)
    if (k % 2 == 0) {
      for (uword l=0; l<4; ++l) {

        const arma::mat &v_mat_ref = value.slice(l);
        arma::mat &v_mat = value_prime.slice(l);

        for (uword i=0; i<par->num_grid_k_uneven; ++i)
          v_mat.row(i) = arma::polyval(arma::polyfit(par->kbar_uneven, v_mat_ref.row(i), par->num_grid_kbar_uneven-1), kbar_prime);
      }
    }

    state->trans = par->pi.row(k);

    for (uword j=0; j<par->num_grid_kbar_even; ++j) {

      // store 4 cubic splines to Params
      for (uword l=0; l<4; ++l) {
        arma::mat &v_mat = value_prime.slice(l);
        state->css[l] = cubic_spline(v_mat.col(j).begin(), v_mat.col(j).end(), par->num_grid_k_uneven);
      }

      for (uword i=0; i<par->num_grid_k_even; ++i) {
        // income
        state->y = par->y_policy(i, j, k);
        inp[0] = 0.;
        success = optim::sumt(inp, val_fn, state, constr_fn, state);
        
        policy(i, j, k) = success ? inp[0] : bisection(state);
      }
    }
  }
}

void policy_to_value(const Params *par,
                     const arma::cube &policy,
                     const arma::cube &value_c,
                     arma::cube &value) {
  // ...
  double discount;
  arma::mat k_prime, c, value_prime(par->num_grid_k_uneven, par->num_grid_kbar_uneven);
  for (uword k=0; k<4; ++k) {
    arma::mat &v_mat = value.slice(k);
    // http://arma.sourceforge.net/docs.html#interp2
    arma::interp2(par->kbar_even, par->k_even, value_c.slice(k), par->kbar_uneven, par->k_uneven, k_prime);

    // compute utility part
    c = par->y_value.slice(k) - k_prime;
    v_mat = utility(c, par->sigma);

    for (uword l=0; l<4; ++k) {
      discount = par->beta * par->pi(k, l);
      const arma::mat &v_mat_ref = value_c.slice(l);
      const arma::rowvec &kbar_prime = par->kbar_prime_value.row(l/2);

      // interpolate at kbar direction (rowwise)
      // http://arma.sourceforge.net/docs.html#polyfit
      for (uword i=0; i<par->num_grid_k_uneven; ++i) 
        value_prime.row(i) = arma::polyval(arma::polyfit(par->kbar_uneven, v_mat_ref.row(i), par->num_grid_kbar_uneven-1), kbar_prime);

      // add discounted future value
      for (uword j=0; j<par->num_grid_kbar_uneven; ++j) {
        cubic_spline cs(value_prime.col(j).begin(), value_prime.col(j).end(), par->num_grid_k_uneven);
        v_mat.col(j) += discount * cs(k_prime.col(j));
      }
    }
  }
}

arma::cube inner_loop(Params *par) {
  // ...
  // value and policy
  arma::cube value(par->num_grid_k_uneven, par->num_grid_kbar_uneven, 4, arma::fill::randu);
  arma::cube policy(par->num_grid_k_even, par->num_grid_kbar_even, 4, arma::fill::randu);

  uword count = 0;
  double err = 1.0;
  arma::cube value_c, policy_c;

  while (err>1e-6 && count<4000) {
    value_c = value;
    policy_c = policy;

    // update policy given value
    // keep policy valid in this function
    value_to_policy(par, value, policy);

    // update value given policy
    policy_to_value(par, policy, value_c, value);

    if (count % 10 == 0)
      err = arma::abs(policy_c - policy).max();

    ++count;
  }

  return policy;
}

// [[Rcpp::export]]
List krusell_smith(double beta,
                   double delta,
                   double sigma,
                   double alpha) {
  /* arguments
  beta : discount factor
  delta: depreciation rate
  sigma: risk aversion
  alpha: capital share
  */
  Params *par = new Params(beta, delta, sigma, alpha);
  Rcout << "11" << std::endl;
  uword count = 0;
  double err = 1.0;
  arma::cube policy;
  while (err>1e-4 && count<80) {
    policy = inner_loop(par);
    break;
  }
  // arma::vec demo(params->num_grid_k_even);

  return List::create(
    _["law_of_motion"] = par->H,
    _["trans_mat"] = par->pi,
    _["policy"] = policy
  );
}

int main() {
  List res = krusell_smith(0.99, 0.025, 1.0, 0.36);
  return 0;
}