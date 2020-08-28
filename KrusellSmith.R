library(Rcpp)

Sys.setenv("PKG_CXXFLAGS"="-I/opt/optimlib/include")
sourceCpp("KrusellSmith.cpp", verbose = TRUE)

res <- krusell_smith(0.99, 0.025, 1, 0.36)

res$law_of_motion

mat <- res$trans_mat

sum(mat[4,])