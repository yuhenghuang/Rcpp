library(Rcpp)

sourceCpp("KalmanFilter.cpp")

pos = matrix(c(1:1000/500, sin(1:1000/500)), ncol = 2, byrow = TRUE)

res <- KalmanFilter(pos)
