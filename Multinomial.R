library(Rcpp)

sourceCpp("Multinomial.cpp")


prob <- runif(200)
prob <- prob/sum(prob) # standardize the probabilities
size <- 500
n <- 20

set.seed(10)
sim_r <- rmultinom(n, size, prob)
set.seed(10)
sim_rcpp <- rmultinom_arma(n, size, prob)
all.equal(sim_r, sim_rcpp)

cores <- parallel::detectCores(logical = FALSE)

microbenchmark::microbenchmark(
  rmultinom(1000, size, prob),
  rmultinom_arma(1000, size, prob),
  rmultinom_omp(1000, size, prob, cores)
)


lambda <- runif(200, 0.5, 3)
set.seed(10)
pois_sim_r <- rpois(length(lambda) + 5, lambda)
set.seed(10)
pois_sim_rcpp <- rpois_arma(length(lambda) + 5, lambda)
all.equal(pois_sim_r, pois_sim_rcpp)
