library(grf)

test_that("quantile forests have reasonable split frequencies", {
  p <- 10
  n <- 500
  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + 100 * (X[, i] > 0))

  qrf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry = p, min.node.size = 10, sample.fraction = 0.632)
  split.frequencies <- split_frequencies(qrf, 4)
  expect_gt(split.frequencies[1, i] / sum(split.frequencies[1, ]), 1 / 2)
})

test_that("quantile forests with regression splitting are identical to regression forests", {
  p <- 10
  n <- 500
  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + 100 * (X[, i] > 0))

  set.seed(1234)
  qrf <- quantile_forest(X, Y,
    quantiles = c(0.1, 0.5, 0.9), regression.splitting = TRUE,
    mtry = p, min.node.size = 10, sample.fraction = 0.632
  )

  set.seed(1234)
  rrf <- regression_forest(X, Y, mtry = p, min.node.size = 10, sample.fraction = 0.632, ci.group.size = 1)

  qrf.split.frequencies <- split_frequencies(qrf, 4)
  rrf.split.frequencies <- split_frequencies(rrf, 4)
  expect_equal(qrf.split.frequencies, rrf.split.frequencies)
})

test_that("quantile forest predictions are positive given positive outcomes", {
  p <- 10
  n <- 500
  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  X.new <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- runif(n) + 100 * (X[, i] > 0)

  qrf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry = p, min.node.size = 10, sample.fraction = 0.632)
  expect_true(all(predict(qrf)$predictions > 0))
  expect_true(all(predict(qrf, X.new)$predictions > 0))
})

test_that("quantile forest predictions for 90th percentile are strongly positively correlated
                     with covariate strongly positively correlated with Y", {
  p <- 10
  n <- 500
  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  X.new <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- runif(n) + 100 * (X[, i] > 0)

  qrf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry = p, min.node.size = 10, sample.fraction = 0.632)
  expect_gt(cor(predict(qrf, quantiles = .9)$predictions, X[, i]), 0.5)
  expect_gt(cor(predict(qrf, X.new, quantiles = .9)$predictions, X.new[, i]), 0.5)
})

test_that("quantile_forest works as expected with missing values", {
  n <- 500
  p <- 5
  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + 100 * (X[, i] > 0))
  quantiles <- c(0.1, 0.5, 0.9)

  nmissing <- 250
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  # MIA with data duplication
  Xl <- X
  Xr <- X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)

  rf.mia <- quantile_forest(X.mia, Y, quantiles = quantiles, seed = 123)
  rf <- quantile_forest(X, Y, quantiles = quantiles, seed = 123)

  mean.diff.oob <- colMeans((predict(rf)$predictions - predict(rf.mia)$predictions))
  mean.diff <- colMeans((predict(rf, X)$predictions - predict(rf.mia, X.mia)$predictions))

  expect_equal(mean.diff.oob, c(0, 0, 0), tolerance = 0.5)
  expect_equal(mean.diff, c(0, 0, 0), tolerance = 0.5)
})
