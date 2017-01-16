rm(list = ls())

setwd("~/git_local/split-relabel/experiments/instrumental_stress_test")

library(gradient.forest)
source("../baselines.R")

n = 1000
n.test = 100
p = 10

X = matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, 1 / (1 + exp(X[,2] / 6)))
Z = W
Y = (2 * W - 1) / 2 * X[,1] + 3 * X[,2] + rnorm(n)

X.test = matrix(rnorm(n.test * p), n.test, p)
tau.true = X.test[,1]

forest.1 = instrumental.forest(X, Y, W, Z, precompute.nuisance = TRUE)
tau.1 = predict(forest.1, newdata = X.test)

forest.2 = instrumental.forest(X, Y, W, Z, precompute.nuisance = FALSE)
tau.2 = predict(forest.2, newdata = X.test)

tau.3 = iv.series(X, Y, W, Z, X.test, interact = FALSE)
tau.4 = iv.series(X, Y, W, Z, X.test, interact = TRUE)
tau.5 = iv.knn(X, Y, W, Z, X.test, k = 100)

c(mean((tau.1 - tau.true)^2), mean((tau.2 - tau.true)^2), mean((tau.3 - tau.true)^2), mean((tau.4 - tau.true)^2), mean((tau.5 - tau.true)^2))

