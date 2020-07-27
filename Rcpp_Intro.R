library(Rcpp)
library(Zelig)

# optional
# redirect cache directory of rcpp
options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")


sourceCpp("./Rcpp_Intro.cpp")

convolveC(1:3, 1:8)


set.seed(883)

initdata <- rnorm(1000, 21, 10)

result_r <- bootstrapC(initdata)


data("turnout")

mY <- matrix(turnout$vote)
mX <- cbind(1,
            turnout$income,
            turnout$educate,
            turnout$age
)


lm(vote~income + educate + age, data=turnout)

fastLM(mX, mY)$coefficient


library(microbenchmark)

microbenchmark(
  seq = (
    lm(vote~income + educate + age, data=turnout)
  ),
  par = (
    fastLM(mX, mY)
  ),
  times = 20
)
