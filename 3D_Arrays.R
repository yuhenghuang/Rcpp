library(Rcpp)

# optional
# redirect cache directory of rcpp
# options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")
# Sys.setenv("PKG_CXXSTD"="-std=c++17")
# Sys.setenv("PKG_CXXFLAGS"="-O3 -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include")
# Sys.setenv("PKG_LIBS"="-L/opt/homebrew/opt/libomp/lib -lomp")

sourceCpp("3D_Arrays.cpp", verbose = TRUE, rebuild = TRUE)


# Matrix Dimensions
xdim <- 200
ydim <- 200

# Number of Time Steps
tdim <- 50

# Generation of 3D Arrays

# Initial Neighborhoods
a <- array(0:1, dim = c(xdim, ydim, tdim))

res_omp <- cube_omp(a, 10)

cube_r_cache <- function(a, xdim, ydim, tdim){
  res = a
  for (t in 1:tdim) {
    temp_time <- a[,,t]
    for (x in 3:(xdim-2)) {
      temp_row <- temp_time[(x-2):(x+2),]
      for (y in 3:(ydim-2)) {
        res[x,y,t] <- sum(temp_row[,(y-2):(y+2)])
      }
    }
  }
  res
}


res_r <- cube_r_cache(a, xdim, ydim, tdim)

all.equal(res_omp, res_r)


res_conv <- cube_conv(a, 10)

all.equal(res_conv, res_r)


microbenchmark::microbenchmark(
  cube_omp(a, 10),
  cube_conv(a, 10),
  cube_r_cache(a, xdim, ydim, tdim)
)
