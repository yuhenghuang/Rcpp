#include <Rcpp.h>
#include <string>
#include <iostream>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
void test_r_func(Function f) {
  String s("r string");
  std::string ss("cpp string");

  f(s);
  f(ss);
}


// [[Rcpp::export]]
void test_cin(Function f) {
  std::string line;
  // the same as stdin in R session
  getline(std::cin, line);

  f(line);
}
