library(Rcpp)

# optional
# redirect cache directory of rcpp
options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

sourceCpp("./Multivariate_Normal_Density.cpp")

cores <- parallel::detectCores(logical = FALSE)

set.seed(123)
sigma <- bayesm::rwishart(10, diag(8))$IW
mu <- rnorm(8)
X <- mvtnorm::rmvnorm(900000, mu, sigma)

library(rbenchmark)

benchmark(
  dmvnorm_arma = dmvnorm_arma(X, mu, sigma),
  dmvnorm_arma_mc = dmvnorm_arma_mc(X, mu, sigma, cores = cores),
  dmvnorm_arma_fast = dmvnorm_arma_fast(X, mu, sigma),
  order = 'relative', replications = 100
)
