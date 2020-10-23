library(RcppMLPACK)

Rcpp::sourceCpp("MLPACK2.cpp")


logistic_regression(matrix(c(1, 2, 3, 1, 2, 3), nrow=2, byrow=TRUE), matrix(c(1L, 1L, 0L), nrow = 1))


data(trainSet)
trainmat <- t(trainSet[, -5]) ## train data
trainlab <- trainSet[, 5] 

linear_regression(trainmat, trainlab)

naive_bayes_classifier(trainmat, trainlab, 2L)


testmat <- t(testSet[, -5])   ## test data
testlab <- testSet[, 5]             
res <- naive_bayes_classifier(trainmat, trainlab, 2L, testmat)   ## also classify
res


all.equal(res[[4]][1,], testlab)


res <- nn_ml(t(randu), 2L)

res$clusters

res$result


res <- kmeans(t(randu), 20L)

dim(res$centroids)

length(res$assignments)
