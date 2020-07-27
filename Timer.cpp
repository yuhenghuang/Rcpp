#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector useTimer() {
  const int n = 1000000;
  Timer timer;
  timer.step("start");

  for (int i=0; i<n; ++i) {
    GetRNGstate();
    PutRNGstate();
  }

  timer.step("get/put");

  for (int i=0; i<n; ++i) {
    GetRNGstate();
    rnorm(10, 0., 1.);
    PutRNGstate();
  }

  timer.step("g/p+rnorm()");
  
  for (int i=0; i<n; ++i) {

  }

  timer.step("empty loop");

  NumericVector res(timer);
  for (int i=0; i<res.size(); ++i)
    res[i] /= n;
    
  return res;
}