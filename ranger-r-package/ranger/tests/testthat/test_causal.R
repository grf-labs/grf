library(ranger)
context("ranger_causal")

p = 10
n = 100
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
W = rbinom(n, 1, 1/(1 + exp(X[,3])))
Y = round(100 * (2*W - 1) * (X[,1] > 0) + X[,2] + rnorm(n), 2)
D = data.frame(X=X, Y=Y, W=W)

test_that("OOB predictions are consistent across training and testing", {
  rf <- ranger(Y ~ ., D, status.variable.name = "W", causal = TRUE, num.trees = 1, write.forest = TRUE)
  oob_samples <- rf$predictions != 0.0
  pred <- predict(rf, D[oob_samples,])
  
  expect_equal(rf$predictions[oob_samples], pred$predictions)
})