library(Rcpp)
library(tidyverse)
library(Zelig)

options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")

data(turnout)

head(turnout)

fit0 <- glm(vote ~ income + educate + age,
            data = turnout,
            family = binomial(link = "probit"))

fit0

sourceCpp("./EM_Algorithm.cpp")

mY <- matrix(turnout$vote)
mX <- cbind(1,
            turnout$income,
            turnout$educate,
            turnout$age
)

fit1 <- em1(
  y = mY,
  X = mX,
  maxit = 20
  )

fit1$beta

head(fit1$eystar)


sourceCpp("./EM_Algorithm.cpp")

fit2 <- em2(
  y = mY,
  X = mX,
  maxit = 20
)

fit2$beta

head(fit2$eystar)


sourceCpp("./EM_Algorithm.cpp")

fit3 <- em3(
  y = mY,
  X = mX,
  maxit = 100
)

fit3$beta

head(fit3$eystar)
