library(Rcpp)

# optional
# redirect cache directory of rcpp
options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp("C++11.cpp")

useAuto()

ls = useInitLists()
ls
typeof(ls)

simpleProd(1:5)
