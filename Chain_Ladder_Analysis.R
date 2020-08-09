library(Rcpp)

sourceCpp("Chain_Ladder_Analysis.cpp")


# Age-To-Age Factors
ageFact <- seq(1.9, 1, by = -.1)

# Inflation Rate
infRate <- 1.02

# Function to reverse matrix columns
revCols <- function(x) {
  x[,ncol(x):1]
}

# Similar to jitter()
shake <- function(vec, sigmaScale = 100) {
  rnorm(n = length(vec), mean = vec, sd = vec/sigmaScale)
}

# Row generation funtion
GenerateRow <- function(iDev, dFactors = cumprod(ageFact), 
                        dInflationRate = 1.02, initClaim = 154) {
  shake(initClaim)*shake(c(1, dFactors))*(dInflationRate^iDev)
}

# Function to generate a claims matrix
GenerateTriangle <- function(iSize, ...) {
  indices = 1:iSize
  mClaimTri = t(sapply(indices, GenerateRow, ...))
  # Reverse columns to get the claims triangle
  mClaimTri = revCols(mClaimTri)
  # Assign nan to lower triangle
  mClaimTri[lower.tri(mClaimTri)] = NA
  mClaimTri = revCols(mClaimTri)
  return(mClaimTri)
}


x <- GenerateTriangle(11)

res_cpp <- getChainSquare_cpp(x)


# Get claims factor at a particular column index
GetFactorR <- function(index, mTri) {
  fact = matrix(mTri[-c((nrow(mTri) - index + 1):nrow(mTri)), index:(index + 1)], ncol = 2)
  fact = c(sum(fact[,1]), sum(fact[,2]))  
  return(fact[2]/fact[1])
}

# Function to carry out Chain Ladder on a claims triangle
GetChainSquareR <- function(mClaimTri) {
  nCols <- ncol(mClaimTri)
  dFactors = sapply(1:(nCols - 1), GetFactorR, mTri = mClaimTri)
  dAntiDiag = diag(revCols(mClaimTri))[2:nCols]
  for(index in 1:length(dAntiDiag)) {
    mClaimTri[index + 1, (nCols - index + 1):nCols] = 
      dAntiDiag[index]*cumprod(dFactors[(nCols - index):(nCols - 1)])
  }
  mClaimTri
}


res_r <- GetChainSquareR(x)

all.equal(res_r, res_cpp)

microbenchmark::microbenchmark(
  GetChainSquareR(x), 
  getChainSquare_cpp(x), 
  times = 10000L)
