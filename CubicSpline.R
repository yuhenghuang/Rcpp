library(Rcpp)

sourceCpp("CubicSpline.cpp")

x <- seq(1, 10 ,by = 0.25)

y <- sin(x)

z <- c(2, 3, 4)

res <- cubic_spline_test(x, y, z)

all.equal(res[,1], sin(z))
