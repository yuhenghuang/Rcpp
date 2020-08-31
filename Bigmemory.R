library(bigmemory)
library(microbenchmark)

Rcpp::sourceCpp("Bigmemory.cpp")

set.seed(4)
ridx <- sample(1:1000, 100)
m <- matrix(rnorm(60*1000), nrow=60)
bigm <- as.big.matrix(m)

options(width=100) # Make sure output looks ok


microbenchmark(res1 <- colSums(m[,ridx]), 
               res2 <- big_col_sum(bigm@address, ridx))

all.equal(res1, res2[1,])
