library(grf)

test_that("instrumental forests give reasonable estimates", {
  p <- 6
  n <- 2000
  n.test <- 2000

  X <- matrix(rnorm(n * p), n, p)

  eps <- rnorm(n)
  Z <- rbinom(n, 1, 2 / 3)
  filter <- rbinom(n, 1, 1 / (1 + exp(-1 * eps)))
  W <- Z * filter

  tau <- apply(X[, 1:2], 1, function(xx) sum(pmax(0, xx)))
  mu <- apply(X[, 2 + 1:2], 1, function(xx) sum(pmax(0, xx)))

  Y <- (2 * W - 1) / 2 * tau + mu + eps

  X.test <- matrix(rnorm(n.test * p), n.test, p)
  tau.true <- apply(X.test[, 1:2], 1, function(xx) sum(pmax(0, xx)))

  forest.iv <- instrumental_forest(X, Y, W, Z, num.trees = 4000)

  preds.iv.oob <- predict(forest.iv, estimate.variance = TRUE)
  preds.iv <- predict(forest.iv, X.test, estimate.variance = TRUE)

  expect_true(all(preds.iv$variance.estimate > 0))
  expect_true(all(preds.iv.oob$variance.estimate > 0))

  error <- preds.iv$predictions - tau.true
  expect_lt(mean(error^2), 0.4)

  error.oob <- preds.iv.oob$predictions - tau
  expect_lt(mean(error.oob^2), 0.4)

  Z <- error / sqrt(preds.iv$variance.estimate)
  Z.oob <- error.oob / sqrt(preds.iv.oob$variance.estimate)

  expect_lt(mean(abs(Z) > 1), 0.75)
  expect_lt(mean(abs(Z.oob) > 1), 0.75)
})

test_that("sample weighted instrumental forest is estimated with kernel weights `forest.weights * sample.weights`", {
  p <- 4
  n <- 500
  X <- matrix(rnorm(n * p), n, p)
  eps <- rnorm(n)
  Z <- rbinom(n, 1, 2 / 3)
  filter <- rbinom(n, 1, 1 / (1 + exp(-1 * eps)))
  W <- Z * filter
  tau <- apply(X[, 1:2], 1, function(xx) sum(pmax(0, xx)))
  mu <- apply(X[, 2 + 1:2], 1, function(xx) sum(pmax(0, xx)))
  Y <- (2 * W - 1) / 2 * tau + mu + eps
  e.cc <- 1 / (1 + exp(-1 * X[, 1]))
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc
  ivf <- instrumental_forest(X, Y, W, Z, Y.hat = 0, W.hat = 0, Z.hat = 0, sample.weights = sample.weights, num.trees = 250)

  x1 <- X[1, , drop = F]
  theta1 <- predict(ivf, x1)$predictions
  alpha1 <- get_forest_weights(ivf, x1)[1, ]
  R <- predict(lm(W ~ Z, weights = alpha1 * sample.weights))
  theta1.lm <- lm(Y ~ R, weights = alpha1 * sample.weights)

  expect_equal(theta1, theta1.lm$coefficients[[2]], tolerance = 1e-10)
})

test_that("instrumental forest predictions and variance estimates are invariant to scaling of the sample weights.", {
  p <- 4
  n <- 200
  X <- matrix(rnorm(n * p), n, p)
  eps <- rnorm(n)
  Z <- rbinom(n, 1, 2 / 3)
  filter <- rbinom(n, 1, 1 / (1 + exp(-1 * eps)))
  W <- Z * filter
  tau <- apply(X[, 1:2], 1, function(xx) sum(pmax(0, xx)))
  mu <- apply(X[, 2 + 1:2], 1, function(xx) sum(pmax(0, xx)))
  Y <- (2 * W - 1) / 2 * tau + mu + eps
  e.cc <- 1 / (1 + exp(-1 * X[, 1]))
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc

  # The multiple is a power of 2 to avoid rounding errors allowing for exact comparison
  # between two forest with the same seed.
  forest.1 <- instrumental_forest(X, Y, W, Z, sample.weights = sample.weights, num.trees = 250, seed = 1)
  forest.2 <- instrumental_forest(X, Y, W, Z, sample.weights = 64 * sample.weights, num.trees = 250, seed = 1)
  pred.1 <- predict(forest.1, estimate.variance = TRUE)
  pred.2 <- predict(forest.2, estimate.variance = TRUE)

  expect_equal(pred.1$predictions, pred.2$predictions, tolerance = 1e-10)
  expect_equal(pred.1$variance.estimates, pred.2$variance.estimates, tolerance = 1e-10)
  expect_equal(pred.1$debiased.error, pred.2$debiased.error, tolerance = 1e-10)
})

test_that("instrumental forests with censoring and ipcc weights compare to forests given full data as you would expect:
 acceptable error, greater variance, and 60% coverage of 95% CIs", {
  p <- 4
  n <- 1000
  n.test <- 200

  X <- matrix(rnorm(n * p), n, p)

  eps <- rnorm(n)
  Z <- rbinom(n, 1, 2 / 3)
  filter <- rbinom(n, 1, 1 / (1 + exp(-1 * eps)))
  W <- Z * filter

  tau <- apply(X[, 1:2], 1, function(xx) sum(pmax(0, xx)))
  mu <- apply(X[, 2 + 1:2], 1, function(xx) sum(pmax(0, xx)))

  Y <- (2 * W - 1) / 2 * tau + mu + eps

  e.cc <- 1 / (1 + exp(-1 * X[, 1]))
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc

  X.test <- matrix(rnorm(n.test * p), n.test, p)
  tau.true <- apply(X.test[, 1:2], 1, function(xx) sum(pmax(0, xx)))

  forest.iv.full <- forest.iv.ipcc <- instrumental_forest(X, Y, W, Z, num.trees = 4000)
  forest.iv.ipcc <- instrumental_forest(X[cc, ], Y[cc], W[cc], Z[cc],
    sample.weights = sample.weights[cc], num.trees = 4000
  )

  preds.iv.full <- predict(forest.iv.full, X.test, estimate.variance = TRUE)
  preds.iv.ipcc <- predict(forest.iv.ipcc, X.test, estimate.variance = TRUE)

  error.full <- preds.iv.full$predictions - tau.true
  error.ipcc <- preds.iv.ipcc$predictions - tau.true
  expect_lt(mean(error.full^2), 0.4)
  expect_lt(mean(error.ipcc^2), 0.6)

  expect_lte(mean(preds.iv.full$variance.estimate) / mean(preds.iv.ipcc$variance.estimate), 1)

  Z.full <- error.full / sqrt(preds.iv.full$variance.estimate)
  Z.ipcc <- error.ipcc / sqrt(preds.iv.ipcc$variance.estimate)

  # ask for a very weak approximation of coverage.
  coverage.full <- mean(abs(Z.full) <= 2)
  coverage.ipcc <- mean(abs(Z.ipcc) <= 2)
  expect_gte(coverage.full, 0.6)
  expect_gte(coverage.ipcc, 0.6)
})
