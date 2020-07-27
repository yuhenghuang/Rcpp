library(Rcpp)

# optional
# redirect cache directory of rcpp
options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")

sourceCpp("Rout.cpp")

ShowValue(1.23)

mat <- matrix(1:9, 3, 3)
showMatrix(mat)

callPrint(1:3)

callPrint(LETTERS[1:3])

callPrint(mat)

callPrint(globalenv())


useOperatorOnVector(1:10)

useOperatorOnMatrix(mat)
