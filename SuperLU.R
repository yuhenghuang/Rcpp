library(Rcpp)

Sys.setenv("PKG_LIBS"="-lsuperlu")

sourceCpp("SuperLU.cpp")

superLU_demo()
