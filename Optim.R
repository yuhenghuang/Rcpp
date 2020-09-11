library(Rcpp)

Sys.setenv("PKG_CXXFLAGS"="-I/opt/optimlib/include/")

sourceCpp("Optim.cpp")

ackley_test()


logistic_test()
