#ifdef _OPENMP
#include <omp.h>
#endif

#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
double long_computation(int nb) {
  double sum = 0.0;
  for (int i=0; i<nb; ++i)
    for (int j=0; j<nb; ++j)
      sum += R::dlnorm(i+j, 0.0, 1.0, 0);

  return sum + nb;
}

// [[Rcpp::export]]
double long_computation2(int nb) {
  double sum = 0.0;
  Progress p(0, false);
  for (int i=0; i<nb; ++i) {
    if (Progress::check_abort())
      return -1.0;
    for (int j=0; j<nb; ++j)
      sum += R::dlnorm(i+j, 0.0, 1.0, 0);
  }
  return sum + nb;
}

// [[Rcpp::export]]
double long_computation3(int nb) {
  double sum = 0.0;
  Progress p(nb*nb, true);
  for (int i=0; i<nb; ++i) {
    if (Progress::check_abort())
      return -1.0;
    for (int j=0; j<nb; ++j) {
      p.increment();
      sum += R::dlnorm(i+j, 0.0, 1.0, 0);
    }
  }
  return sum + nb;
}

// [[Rcpp::export]]
double long_computation_omp(int nb, int threads=1) {
#ifdef _OPENMP
  if (threads>0)
    omp_set_num_threads(threads);
  Rcout << "Num. of threads = " << omp_get_max_threads() << std::endl;
#endif

  double sum = 0.0, thread_sum;
#pragma omp parallel for reduction(+ : sum) private(thread_sum)
  for (int i=0; i<nb; ++i) {
    thread_sum = 0;
    for (int j=0; j<nb; ++j)
      thread_sum += R::dlnorm(i+j, 0.0, 1.0, 0);
    sum += thread_sum;
  }

  return sum + nb;
}

// [[Rcpp::export]]
double long_computation_omp2(int nb, int threads=1) {
#ifdef _OPENMP
  if (threads>0)
    omp_set_num_threads(threads);
  Rcout << "Num. of threads = " << omp_get_max_threads() << std::endl;
#endif

  Progress p(nb, true);
  double sum = 0.0, thread_sum;
#pragma omp parallel for reduction(+ : sum) private(thread_sum)
  for (int i=0; i<nb; ++i) {
    thread_sum = 0;
    if (!Progress::check_abort()){
      p.increment();
      for (int j=0; j<nb; ++j)
        thread_sum += R::dlnorm(i+j, 0.0, 1.0, 0);
    }
    sum += thread_sum;
  }

  return sum + nb;
}
