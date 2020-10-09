library(grf)

test_that("multi_regression_forest works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)
  nmissing <- 50
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  # A regression forest trained on Y is the same as a multi regression forest
  # trained on [Y 0]
  rf <- regression_forest(X, Y, num.trees = 100, ci.group.size = 1,
                          min.node.size = 1, alpha = 0,
                          seed = 42)
  mrf <- multi_regression_forest(X, cbind(Y, rep(0, n)), num.trees = 100,
                                 min.node.size = 1, alpha = 0,
                                 seed = 42)

  expect_equal(predict(rf)$predictions, predict(mrf)$predictions[, 1])
  expect_equal(predict(rf, X)$predictions, predict(mrf, X)$predictions[, 1])

  expect_equal(ncol((mrf)$predictions), 2)

  # A multi regression forest trained on duplicated outcomes yields the same result
  X <- matrix(rnorm(n * p), n, p)
  YY <- 2 * X[, 1:2] + rnorm(n)
  mrf <- multi_regression_forest(X, YY, num.trees = 100, seed = 1)
  mrf.dup <- multi_regression_forest(X, cbind(YY, YY), num.trees = 100, seed = 1)
  expect_equal(predict(mrf)$predictions, predict(mrf.dup)$predictions[, 1:2])
})

test_that("multi_regression_forest is similar to two regression_forest", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <-  X[, 1, drop = F] %*% cbind(1, 2) + rnorm(n)

  rf1 <- regression_forest(X, Y[, 1], num.trees = 1000)
  rf2 <- regression_forest(X, Y[, 2], num.trees = 1000)
  mrf <- multi_regression_forest(X, Y, num.trees = 1000)

  diff1 <- mean((predict(mrf)$predictions[, 1] - predict(rf1)$predictions)^2)
  diff2 <- mean((predict(mrf)$predictions[, 2] - predict(rf2)$predictions)^2)

  expect_true(diff1 < 0.01)
  expect_true(diff2 < 0.01)
})
