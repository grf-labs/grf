rm(list = ls())

library(randomForest)
library(gradient.forest)
library(Rcpp)

p = 20
n = 1000

X = round(matrix(runif(n * p, -1, 1), n, p), 2)
Y = round(X[,1] + rnorm(n), 2)

X.test = matrix(0, 41, p)
X.test[,1] = seq(-1, 1, by = 0.05)

plot(X.test[,1], predict(randomForest(X, Y, replace=FALSE), X.test))
points(X.test[,1], predict(regression.forest(X, Y, sample.fraction=0.632, num.threads=4), X.test, num.threads=4), col = 2)
abline(0, 1)
