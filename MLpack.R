library(Rcpp)
library(rbenchmark)

sourceCpp("MLpack.cpp")

data(wine, package = "HDclassif")

data(trees, package = "datasets")

set.seed(10)
res_k_cpp <- as.numeric(kmeans_ml(t(trees), 3)$result)
set.seed(10)
res_k <- kmeans(trees, 3)$cluster

all.equal(res_k, res_k_cpp+1)

wine$class <- NULL

mat <- t(wine)

benchmark(
  kmeans_ml(mat, 3),
  kmeans(wine, 3),
  order = "relative"
)

data("trainSet", package = "RcppMLPACK")
mat <- t(trainSet[, -5])
lab <- trainSet[, 5]
logistic_regression(mat, lab)


X <- with(trees, cbind(log(Girth), log(Height)))
y <- with(trees, log(Volume))

mod <- lm(y ~ X)

mlfit <- lr_cpp(X, y)
