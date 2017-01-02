setwd("/Users/jtibshirani/code/split-relabel/r-package")
source("build_package.R")

quantiles = c(0.25, 0.5, 0.75)

p = 20
n = 100

X = matrix(2 * runif(n * p) - 1, n, p)
Y = rnorm(n) * (1 + (X[,1] > 0))
D = data.frame(X=X, Y=Y)

M = as.matrix(sapply(D, as.numeric))

forest <- quantile_train(quantiles, M, 20, matrix(NA, 0, 0), character(0), 5, 4, FALSE, 2, 5, TRUE, FALSE, 0.63, numeric(0), 42)

quantile_predict(forest, quantiles, M, matrix(NA, 0, 0), character(0), 3)