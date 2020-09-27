library(Rcpp)

sourceCpp("Function.cpp")

fn <- function(s) {
  print(typeof(s))
  print(s)
}

test_r_func(fn)

# test_cin(fn)

# x <- readLines("stdin", 1)
