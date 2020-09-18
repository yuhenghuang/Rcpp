library(Rcpp)

sourceCpp("Heap.cpp")


x <- rnorm(1000)

res_cpp <- top_index(x, 30L)

res_r <- tail(order(x), 30L)

all.equal(res_cpp[,1], res_r)

top <- function(x, n) {
  tail(order(x), n)
}


library(microbenchmark)

x <- rnorm(1e5)
microbenchmark( 
  R_order = top(x, 100),
  cpp2    = top_index( x, 100 )
)
