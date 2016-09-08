library(ranger)
context("ranger_quantile")

p = 10
n = 100
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
Y = round(100 * (2*W - 1) * (X[,1] > 0) + X[,2] + rnorm(n), 2)
D = data.frame(X=X, Y=Y)

test_that("Quantile predictions have right dimensions", {
  rf <- ranger(Y ~ ., D, quantile = TRUE, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, D)
})