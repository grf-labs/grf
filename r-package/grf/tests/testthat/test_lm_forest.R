test_that("lm_forest with single W ~ causal forest", {
  # These tests are not done with an epsilon tolerance due to forests'
  # discontinuous nature. Even though these two calls are in principle identical (with the same seed),
  # some splits futher down the tree might deviate by chance due to minor numerical differences in
  # implementation, thus leading to final point predictions that can differ more than an epsilon.
  # For tests that locks in an equivalence between Causal Forest and its multivariate extension
  # see `MultiCausalSplittingRuleTest.cpp`

  # Binary W
  n <- 1500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  wts <- sample(1:2, n, TRUE)

  Y.hat <- predict(regression_forest(X, Y, num.trees = 500, sample.weights = wts))$predictions
  cf <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = 0.5, sample.weights = wts, num.trees = 500, stabilize.splits = FALSE)
  lmf <- lm_forest(X, Y, W, Y.hat = Y.hat, W.hat = 0.5, sample.weights = wts, num.trees = 500)
  expect_lt(mean((predict(cf)$predictions - predict(lmf)$predictions[,,])^2), 0.03)
  expect_equal(mean(predict(cf)$predictions), mean(predict(lmf)$predictions[,,]), tolerance = 0.03)

  # Continuous W
  W <- runif(n)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  Y.hat <- predict(regression_forest(X, Y, num.trees = 500, sample.weights = wts))$predictions
  cfw <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = 0.5, sample.weights = wts, num.trees = 500, stabilize.splits = FALSE)
  lmfw <- lm_forest(X, Y, W, Y.hat = Y.hat, W.hat = 0.5, sample.weights = wts, num.trees = 500)
  expect_lt(mean((predict(cfw)$predictions - predict(lmfw)$predictions[,,])^2), 0.05)
  expect_equal(mean(predict(cfw)$predictions), mean(predict(lmfw)$predictions[,,]), tolerance = 0.05)
})

test_that("lm_forest with dummy W = multi arm causal forest", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  Y <- X[, 1] + 1.5 * (W == "A") + 2.8 * (W == "B") - 4 * (W == "C") + rnorm(n)
  wts <- sample(1:2, n, TRUE)

  Y.hat <- predict(multi_regression_forest(X, Y, num.trees = 250))$predictions
  W.hat <- predict(probability_forest(X, W, num.trees = 250))$predictions

  W.matrix <- model.matrix(~ W - 1)

  mcf <- multi_arm_causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.trees = 250, sample.weights = wts, seed = 42, stabilize.splits = FALSE)
  lmf <- lm_forest(X, Y, W.matrix[, -1], Y.hat = Y.hat, W.hat = W.hat[, -1], num.trees = 250, sample.weights = wts, seed = 42)

  expect_equal(unname(predict(lmf)$predictions), unname(predict(mcf)$predictions))
})

test_that("lm_forest gradient.weights option works as expected", {
  n <- 250
  p <- 5
  K <- 2
  X <- matrix(rnorm(n * p), n, p)
  W <- matrix(runif(n * K), n, K)
  Y <- X[, 1] - W[, 1] * pmax(X[, 2], 0) + W[, 2] + rnorm(n)

  lmf <- lm_forest(X, Y, W, num.trees = 250, gradient.weights = c(0.5, 1), seed = 42)
  lmf2 <- lm_forest(X, cbind(Y, Y), W, num.trees = 250, gradient.weights = c(0.5, 1), seed = 42)
  expect_equal(predict(lmf)$predictions[,,], predict(lmf2)$predictions[,, 1])
  expect_equal(predict(lmf)$predictions[,,], predict(lmf2)$predictions[,, 2])

  lmf3 <- lm_forest(X, Y, W, num.trees = 250, seed = 42)
  lmf4 <- lm_forest(X, Y, W, num.trees = 250, gradient.weights = c(0.5, 0.5), seed = 42)
  expect_equal(predict(lmf3)$predictions, predict(lmf4)$predictions)
})
