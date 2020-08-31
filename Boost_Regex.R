library(Rcpp)

Sys.setenv("PKG_LIBS"="-lboost_regex")

sourceCpp("Boost_Regex.cpp")

s <- c("0000111122223333", "0000 1111 2222 3333", "0000-1111-2222-3333", "000-1111-2222-3333")
regex_demo(s)
