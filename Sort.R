library(Rcpp)

sourceCpp("Sort.cpp")

z <- rnorm(100)

y <- stl_sort(z)

stl_sort_inplace(z)

z[1:10]

