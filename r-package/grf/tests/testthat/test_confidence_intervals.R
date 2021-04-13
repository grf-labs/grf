library(grf)

set.seed(123)

test_that("using big ci.group.size doesn't result in tiny CIs", {
  p <- 4
  n <- 200

  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- (X[, 1] > 0) + rnorm(n)
  X.test <- matrix(2 * runif(2000 * p) - 1, 2000, p)

  forest.big.group <- regression_forest(X, Y, num.trees = 2000, ci.group.size = 200)
  pred.big.group <- predict(forest.big.group, X.test, estimate.variance = TRUE)
  var.big.group <- sqrt(pred.big.group$variance.estimates)

  forest.sm.group <- regression_forest(X, Y, num.trees = 2000, ci.group.size = 10)
  pred.sm.group <- predict(forest.sm.group, X.test, estimate.variance = TRUE)
  var.sm.group <- sqrt(pred.sm.group$variance.estimates)

  expect_gt(min(var.big.group), 0.01)
  expect_gt(min(var.sm.group), 0.01)
  expect_gt(mean(var.big.group) / mean(var.sm.group), 0.8)
})


test_that("regression CIs are reasonable", {
  p <- 6
  n <- 1000
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- (X[, 1] > 0) + 2 * rnorm(n)
  forest <- regression_forest(X, Y)
  preds.oob <- predict(forest, estimate.variance = TRUE)
  error.standardized <- (preds.oob$predictions - (X[, 1] > 0)) / sqrt(preds.oob$variance.estimates)
  expect_lte(mean(abs(error.standardized) > qnorm(0.975)), 0.15)
})

test_that("instrumental CIs are reasonable", {
  n <- 1000
  n.test <- 1000
  p <- 4
  alpha.mu <- 3
  alpha.tau <- 1
  k.tau <- 2
  k.mu <- 2

  X <- matrix(rnorm(n * p), n, p)
  eps <- 2 * rnorm(n)
  Z <- rbinom(n, 1, 2 / 3)
  filter <- rbinom(n, 1, 1 / (1 + exp(-eps)))
  W <- Z * filter
  tau <- alpha.tau * apply(X[, 1:k.tau], 1, function(xx) sum(pmax(0, xx)))
  mu <- alpha.mu * apply(X[, k.tau + 1:k.mu], 1, function(xx) sum(pmax(0, xx)))
  Y <- (2 * W - 1) / 2 * tau + mu + eps

  X.test <- matrix(rnorm(n.test * p), n.test, p)
  tau.true <- alpha.tau * apply(X.test[, 1:k.tau], 1, function(xx) sum(pmax(0, xx)))

  forest <- instrumental_forest(X, Y, W, Z)
  tau.hat <- predict(forest, newdata = X.test, estimate.variance = TRUE)
  error.standardized <- (tau.hat$predictions - tau.true) / sqrt(tau.hat$variance.estimates)
  expect_lte(mean(abs(error.standardized) > qnorm(0.975)), 0.18)
})

test_that("instrumental CIs are invariant to scaling Z", {
  n <- 2000
  p <- 5
  n.test <- 2000
  X <- matrix(rnorm(n * p), n, p)
  W <- rnorm(n)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

  X.test <- matrix(rnorm(n * p), n, p)
  tau.true <- pmax(X.test[, 1], 0)

  y.forest <- regression_forest(X, Y)
  Y.hat <- predict(y.forest)$predictions

  forest <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = 0)
  tau.hat <- predict(forest, newdata = X.test, estimate.variance = TRUE)
  error.standardized <- (tau.hat$predictions - tau.true) / sqrt(tau.hat$variance.estimates)
  expect_lte(mean(abs(error.standardized) > qnorm(0.975)), 0.15)
  expect_gte(mean(abs(error.standardized) > qnorm(0.975)), 0.005)

  Z <- 0.00000001 * W
  forest.iv <- instrumental_forest(X, Y, W, Z,
    Y.hat = Y.hat, W.hat = 0, Z.hat = 0
  )
  tau.hat.iv <- predict(forest.iv, newdata = X.test, estimate.variance = TRUE)
  error.standardized.iv <- (tau.hat.iv$predictions - tau.true) / sqrt(tau.hat.iv$variance.estimates)
  expect_lte(mean(abs(error.standardized.iv) > qnorm(0.975)), 0.15)
  expect_gte(mean(abs(error.standardized.iv) > qnorm(0.975)), 0.005)

  kst <- ks.test(error.standardized, error.standardized.iv)
  expect_lte(kst$statistic, 0.05)
})

test_that("LL causal CIs are reasonable", {
   n <- 1000
   n.test <- 1000
   p <- 4
   X <- matrix(rnorm(n * p), n, p)
   W <- rbinom(n, 1, 0.5)
   TAU <- 2 * X[, 1]
   Y <- W * TAU + 0.5 * rnorm(n)

   forest <- causal_forest(X, Y, W)
   tau.hat <- predict(forest, linear.correction.variables = 1:ncol(X), estimate.variance = TRUE)
   error.standardized <- (tau.hat$predictions - TAU) / sqrt(tau.hat$variance.estimates)
   expect_lt(mean(abs(error.standardized) > qnorm(0.975)), 0.1)

   X.test <- matrix(rnorm(n.test * p), n.test, p)
   TAU.test <- 2 * X.test[, 1]

   tau.hat.test <- predict(forest, X.test, linear.correction.variables = 1:ncol(X), estimate.variance = TRUE)
   error.standardized.test <- (tau.hat.test$predictions - TAU.test) / sqrt(tau.hat.test$variance.estimates)
   expect_lt(mean(abs(error.standardized.test) > qnorm(0.975)), 0.1)
})

test_that("LL causal CIs are shorter than standard causal CIS", {
   n <- 1000
   p <- 4
   X <- matrix(rnorm(n * p), n, p)
   W <- rbinom(n, 1, 0.5)
   TAU <- 2 * X[, 1] + X[, 2]
   Y <- W * TAU + rnorm(n)

   forest <- causal_forest(X, Y, W)
   tau.hat.ll <- predict(forest, linear.correction.variables = 1:2, ll.lambda = 0.01, estimate.variance = TRUE)
   tau.hat <- predict(forest, estimate.variance = TRUE)

   preds <- tau.hat$predictions
   vars <- tau.hat$variance.estimates
   lowers <- preds - 1.96 * sqrt(vars)
   uppers <- preds + 1.96 * sqrt(vars)
   coverage.grf <- mean(lowers <= TAU & TAU <= uppers)
   length.grf <- mean(abs(uppers - lowers))

   preds.ll <- tau.hat.ll$predictions
   vars.ll <- tau.hat.ll$variance.estimates
   lowers.ll <- preds.ll - 1.96 * sqrt(vars.ll)
   uppers.ll <- preds.ll + 1.96 * sqrt(vars.ll)
   coverage.ll <- mean(lowers.ll <= TAU & TAU <= uppers.ll)
   length.ll <- mean(abs(uppers.ll - lowers.ll))

   expect_lt(coverage.grf / coverage.ll, 1)
   expect_lt(length.ll / length.grf, 0.8)
})
