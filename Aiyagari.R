library(Rcpp)

sourceCpp("Aiyagari.cpp")


st <- Sys.time()

res <- Aiyagari(0.6, 3.0, 0.4, 4)

total <- Sys.time() - st
print(total)

res$knew

res$k

res$r

res$w
