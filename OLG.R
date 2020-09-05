library(Rcpp)

sourceCpp("OLG.cpp")

res = OLG_model(0.36, 0.1, 1.01, 60)

res$a

res$c

sum(res$c)
