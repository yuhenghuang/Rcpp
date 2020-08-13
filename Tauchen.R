library(Rcpp)

sourceCpp("Tauchen.cpp")

res <- tauchen_test(0.6,
                    0.4*sqrt(1-0.36))

res$grid

res$station

res$pi

sum(res$station)
