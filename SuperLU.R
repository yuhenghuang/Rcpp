library(Rcpp)
library(glue)

Sys.setenv("PKG_LIBS"=glue("{Sys.getenv('PKG_LIBS')} -lsuperlu"))

sourceCpp("SuperLU.cpp", verbose = TRUE)

superLU_demo()
