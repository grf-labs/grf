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
  mrf <- multi_regression_forest(X, cbind(Y, 0), num.trees = 100,
                                 min.node.size = 1, alpha = 0,
                                 seed = 42)

  expect_equal(predict(rf)$predictions, predict(mrf)$predictions[, 1])
  expect_equal(predict(rf, X)$predictions, predict(mrf, X)$predictions[, 1])
  expect_equal(ncol((mrf)$predictions), 2)

  # The above logic holds "symmetrically":

  # A regression forest trained on Y is the same as a multi regression forest
  # trained on [0 Y]
  mrf <- multi_regression_forest(X, cbind(0, Y), num.trees = 100,
                                 min.node.size = 1, alpha = 0,
                                 seed = 42)

  expect_equal(predict(rf)$predictions, predict(mrf)$predictions[, 2])
  expect_equal(predict(rf, X)$predictions, predict(mrf, X)$predictions[, 2])
  expect_equal(ncol((mrf)$predictions), 2)

  # A regression forest trained on Y is the same as a multi regression forest
  # trained on [0 0 Y 0 0 0]
  mrf <- multi_regression_forest(X, cbind(0, 0, Y, 0, 0, 0), num.trees = 100,
                                 min.node.size = 1, alpha = 0,
                                 seed = 42)

  expect_equal(predict(rf)$predictions, predict(mrf)$predictions[, 3])
  expect_equal(predict(rf, X)$predictions, predict(mrf, X)$predictions[, 3])
  expect_equal(ncol((mrf)$predictions), 6)

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

  expect_lt(diff1, 0.01)
  expect_lt(diff2, 0.01)
})

test_that("multi_regression_forest is on parity with regression_forest", {
  # Test that the R package implementation keeps a multi regression forest with one outcome
  # identical to a regression forest.
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <-  X[, 1] * rnorm(n)
  nmissing <- sample(c(1, 150, 400), size = 1)
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  mrf <- multi_regression_forest(X, Y, num.trees = 500, seed = 42)
  rf <- regression_forest(X, Y, num.trees = 500, ci.group.size = 1, seed = 42)

  expect_equal(predict(mrf)$predictions[,], predict(rf)$predictions)
  expect_equal(predict(mrf, X)$predictions[,], predict(rf, X)$predictions)
})

test_that("multi_regression_forest is well calibrated", {
  n <- 250
  p <- 25
  X <- matrix(rnorm(n * p), n, p)

  signal <- pmax(X[, 1], 0)
  nY <- 5
  mu <- matrix(rep(signal, nY), n, nY)
  Y <- mu + matrix(rnorm(n * nY), n, nY)

  # A multi regression forest trained on the same signal with idiosyncratic noise
  # should yield a lower MSE than multiple regression forests
  rf.pred <- apply(Y, 2, function(y) predict(regression_forest(X, y, num.trees = 250))$predictions)
  mrf.pred <- predict(multi_regression_forest(X, Y, num.trees = 250))$predictions

  mse.rf <- mean((rf.pred - mu)^2)
  mse.mrf <- mean((mrf.pred - mu)^2)

  expect_lt(mse.mrf / mse.rf, 0.8)
})

test_that("sample weighted multi_regression_forest is estimated with kernel weights `forest.weights * sample.weights`", {
  n <- 2000
  p <- 5
  obs.prob <- 1 / 20
  Y0 <- rbinom(n, 1, obs.prob / (1 + obs.prob))
  Y <- Y0 + matrix(rnorm(n * 2), n, 2) * 0.01
  X <- matrix(rnorm(n * p), n, p)
  sample.weights <- 1 + Y0 * (1 / obs.prob - 1)
  mrf <- multi_regression_forest(X, Y, sample.weights = sample.weights, num.trees = 250)

  x1 <- X[1, , drop = F]
  theta1 <- predict(mrf, x1)$predictions
  alpha1 <- get_forest_weights(mrf, x1)[1, ]
  theta1.lm <- lm(Y ~ 1, weights = alpha1 * sample.weights)

  expect_equal(unname(theta1[1, ]), theta1.lm$coefficients[1, ], tolerance = 1e-10)
})
