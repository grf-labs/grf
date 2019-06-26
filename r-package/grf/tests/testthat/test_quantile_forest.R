library(grf)

test_that("quantile forests have reasonable split frequencies", {
  p <- 10
  n <- 500
  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + 100 * (X[, i] > 0))

  qrf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry = p, min.node.size = 10, sample.fraction = 0.632)
  split.frequencies <- split_frequencies(qrf, 4)
  expect_true(split.frequencies[1, i] / sum(split.frequencies[1, ]) > 1 / 2)
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
  expect_true(all(predict(qrf) > 0))
  expect_true(all(predict(qrf, X.new) > 0))
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
  expect_true(cor(predict(qrf, quantiles = .9), X[, i]) > .5)
  expect_true(cor(predict(qrf, X.new, quantiles = .9), X.new[, i]) > .5)
})
