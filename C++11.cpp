#include <Rcpp.h>
#include <vector>
#include <string>
// [[Rcpp::plugins(cpp11)]]

using namespace std;

// [[Rcpp::export]]
int useAuto() {
  auto val = 42;
  return val;
}

// [[Rcpp::export]]
vector<string> useInitLists() {
  vector<string> vec = {"larry", "curly", "moe"};
  return vec;
}

// [[Rcpp::export]]
int simpleProd(vector<int> vec) {
  int prod = 1;
  for (int &x : vec) {
    prod *= x;
  }
  return prod;
}