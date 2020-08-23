library(Rcpp)
library(modelr)
library(lmtest)
library(sandwich)
library(plm)

sourceCpp("ClusteredSE.cpp")

data(diamonds, package = "ggplot2")

# X = cbind(1, diamonds$carat)
X <- data.matrix(model_matrix(diamonds, price ~ carat + factor(cut) - 1))

Y = diamonds$price

G = as.integer(diamonds$cut) - 1

res <- cluster_regression(X, Y, G)

res

cluster_raw_r <- function() {
  mod <- plm(price ~ carat + factor(cut) - 1, data = diamonds, index = "cut")
  
  print(coeftest(mod, vcov = vcovHC(mod, cluster = "group")))
  
  return(mod)
}

cluster_raw_r()

library(rbenchmark)

benchmark(
  rcpp = cluster_regression(X, Y, G),
  # rraw = cluster_raw_r(),
  order = "relative", replications = 100
)