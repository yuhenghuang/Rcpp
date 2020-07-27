library(Rcpp)

# optional
# redirect cache directory of rcpp
options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")

sourceCpp("Timer.cpp")

res <- useTimer()
res

diff(res)
