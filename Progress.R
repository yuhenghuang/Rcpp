library(Rcpp)

sourceCpp("Progress.cpp")

system.time(s  <- long_computation(1000))

system.time(s2  <- long_computation2(5000)) # interrupt me

system.time(s3  <- long_computation3(5000))

system.time(s4 <- long_computation_omp(5000, 4))

system.time(s4 <- long_computation_omp2(10000, 4))
