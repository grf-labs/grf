library(grf)

set.seed(3141)

test_that("causal forests can split on the last parameter", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * (X[, 1] + X[, 6]) + rnorm(n)

  forest <- causal_forest(X, Y, W, compute.oob.predictions = FALSE)
  split.freq <- split_frequencies(forest, 10)

  expect_gt(sum(split.freq[, 6]), 0)
})

test_that("causal forests have reasonable split frequencies", {
  n <- 100
  p <- 7
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.2)
  Y <- 1000 * (X[, p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase imbalance.penalty to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc <- causal_forest(X, Y, W, mtry = p, imbalance.penalty = 1.0, stabilize.splits = TRUE, min.node.size = 2)
  split.freq <- split_frequencies(ccc, 4)
  expect_gt(split.freq[1, p] / sum(split.freq[1, ]), 2 / 3)
})

test_that("causal forests without stable splitting have reasonable split frequencies", {
  n <- 100
  p <- 7
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.2)
  Y <- 1000 * (X[, p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase imbalance.penalty to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc <- causal_forest(X, Y, W, mtry = p, imbalance.penalty = 1.0, stabilize.splits = FALSE, min.node.size = 2)
  split.freq <- split_frequencies(ccc, 4)
  expect_gt(split.freq[1, p] / sum(split.freq[1, ]), 2 / 3)
})

test_that("causal forests with a positive imbalance.penalty have reasonable tree depths", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(p = 0.5, size = 1, n = n)
  Y <- 0.5 * X[, 1] * (2 * W - 1) + 0.1 * rnorm(n)

  forest <- causal_forest(X, Y, W, imbalance.penalty = 0.001)
  split.freq <- split_frequencies(forest)
  expect_gt(sum(split.freq[3, ]), 0)
})

test_that("causal forests with a very small imbalance.penalty behave similarly to unpenalized forests.", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- X[, 1] * (2 * W - 1) + 0.1 * rnorm(n)

  forest <- causal_forest(X, Y, W, imbalance.penalty = 0.0)
  forest.large.penalty <- causal_forest(X, Y, W, imbalance.penalty = 100.0)
  forest.small.penalty <- causal_forest(X, Y, W, imbalance.penalty = 1e-6)

  diff.large.penalty <- abs(forest.large.penalty$debiased.error - forest$debiased.error)
  diff.small.penalty <- abs(forest.small.penalty$debiased.error - forest$debiased.error)
  expect_lt(mean(diff.small.penalty), 0.10 * mean(diff.large.penalty))
})

test_that("causal forests behave reasonably with a low treatment probability", {
  n <- 1000
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.1)
  tau <- 0.1
  Y <- X[, 1] + X[, 2] + tau * W + rnorm(n)

  forest <- causal_forest(X, Y, W, stabilize.splits = TRUE)
  tau.hat <- predict(forest)$predictions
  expect_lt(sqrt(mean((tau.hat - tau)^2)), 0.20)
})

test_that("causal forests behave reasonably with small sample size", {
  p <- 5
  n <- 50
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- 100 * (X[, 1] > 0)
  Y <- tau * (W - 0.5) + rnorm(n)
  forest <- causal_forest(X, Y, W,
    stabilize.splits = TRUE,
    min.node.size = 1, mtry = p,
    compute.oob.predictions = FALSE
  )

  X.test <- matrix(rnorm(n * p), n, p)
  tau.test <- 100 * (X.test[, 1] > 0)
  tau.hat <- predict(forest, X.test)$predictions
  expect_lt(sqrt(mean((tau.hat - tau.test)^2)) / 100, 1 / 3)
})

test_that("local linear causal forests work in a simple case", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  TAU <- 2 * X[, 1] + X[, 2]
  Y <- W * TAU + rnorm(n)

  forest <- causal_forest(X, Y, W, num.trees = 400)
  preds.ll <- predict(forest, X, linear.correction.variables = 1:2, ll.lambda = 0.01)
  error.ll <- mean((preds.ll$predictions - TAU)^2)

  preds.rf <- predict(forest, X)
  error.rf <- mean((preds.rf$predictions - TAU)^2)

  expect_lt(error.ll, 0.8 * error.rf)
})

test_that("local linear causal forests with large lambda are equivalent to causal forests", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- 2 * X[, 1] + X[, 2]
  Y <- W * tau + rnorm(n)

  forest <- causal_forest(X, Y, W, num.trees = 400)
  preds.ll <- predict(forest, X, linear.correction.variables = 1:ncol(X), ll.lambda = 1e5)$predictions
  preds.cf <- predict(forest)$predictions

  expect_lt(mean((preds.ll - preds.cf)^2), 0.02)
})

test_that("predictions are invariant to scaling of the sample weights.", {
  n <- 100
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-2 * X[, 1] + 3 * X[, 2]))
  W <- rbinom(n, 1, e)
  tau <- 2 * X[, 1] + X[, 2]
  Y <- W * tau + rnorm(n)
  e.cc <- 1 / (1 + exp(-2 * X[, 1]))
  sample.weights <- 1 / e.cc

  forest.1 <- causal_forest(X, Y, W, sample.weights = sample.weights, seed = 1)
  # The multiple is a power of 2 to avoid rounding errors.
  forest.2 <- causal_forest(X, Y, W, sample.weights = 64 * sample.weights, seed = 1)
  expect_true(all(abs(forest.2$predictions - forest.1$predictions) < 1e-10))
})

test_that("IPCC weighting in the training of a causal forest with missing data improves its complete-data MSE.", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-X[, 1] + 2 * X[, 2]))
  W <- rbinom(n, 1, e)
  tau <- 2 * (X[, 1] > 0 & X[, 5] > 0) -
      0.5 * (X[, 2] > 0) - 0.5 * (X[, 3] > 0) - 0.5 * (X[, 4] > 0)
  Y <- W * tau + rnorm(n)

  e.cc <- 1 - 0.9 * (X[, 1] > 0 & X[, 5] > 0)
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc

  num.trees <- 500
  mse <- function(f) {
      tau.hat <- rep(NA, n)
      tau.hat[cc] <- predict(f)$predictions
      tau.hat[!cc] <- predict(f, X[!cc, ])$predictions
      mean((tau.hat - tau)^2)
  }

  forest <- causal_forest(X[cc, ], Y[cc], W[cc], num.trees = num.trees)
  weighted.forest <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights = sample.weights[cc], num.trees = num.trees)
  expect_lt(mse(weighted.forest) / mse(forest), .9)

  boosted.forest <- causal_forest(X[cc, ], Y[cc], W[cc], orthog.boosting = TRUE, num.trees = num.trees)
  boosted.weighted.forest <- causal_forest(
      X[cc, ], Y[cc], W[cc],
      sample.weights = sample.weights[cc],
      orthog.boosting = TRUE, num.trees = num.trees
  )
  expect_lt(mse(boosted.weighted.forest) / mse(boosted.forest), .9)
})

test_that("Weighting is roughly equivalent to replication of samples", {
  n <- 500
  p <- 2
  num.trees <- 1000

  X <- matrix(rnorm(n * p), n, p)
  e.a <- 1 / (1 + exp(+1 * X[, 1] + 1 * X[, 2]))
  e.b <- 1 / (1 + exp(-1 * X[, 1] - 1 * X[, 2]))
  tau.a <- +1000 + 2 * X[, 1] + X[, 2]
  tau.b <- -1000 + X[, 1] - X[, 2]
  W.a <- rbinom(n, 1, e.a)
  W.b <- rbinom(n, 1, e.b)
  Y.a <- W.a * tau.a + rnorm(n)
  Y.b <- W.b * tau.b + rnorm(n)

  X.rep <- X[rep(1:(n / 2), 5), ]
  W.rep <- c(W.a[rep(1:(n / 2), 4)], W.b[1:(n / 2)])
  Y.rep <- c(Y.a[rep(1:(n / 2), 4)], Y.b[1:(n / 2)])
  sample.weights <- c(rep(4, n / 2), rep(1, n / 2))
  train <- 1:(n / 2)
  test <- (n / 2 + 1):n

  ## regression forest for propensity score
  regression.forest.rep <- regression_forest(X.rep, W.rep, num.trees = num.trees)
  regression.forest.weight <- regression_forest(X[rep(train, 2), ], c(W.a[train], W.b[train]),
    sample.weights = sample.weights, num.trees = num.trees
  )
  regression.forest.biased <- regression_forest(X[rep(train, 2), ], c(W.a[train], W.b[train]),
    num.trees = num.trees
  )

  z_scores <- function(a, b) {
    # conservative. use a variance upper bound
    abs(a$predictions - b$predictions) / sqrt(2 * (a$variance + b$variance))
  }

  expect_true(mean(z_scores(
    predict(regression.forest.rep, X[test, ], estimate.variance = TRUE),
    predict(regression.forest.weight, X[test, ], estimate.variance = TRUE)
  ) <= 1) >= .5)
  expect_true(mean(z_scores(
    predict(regression.forest.rep, X[test, ], estimate.variance = TRUE),
    predict(regression.forest.biased, X[test, ], estimate.variance = TRUE)
  ) >= 1) >= .5)
  ## causal forest
  causal.forest.rep <- causal_forest(X.rep, Y.rep, W.rep,
    W.hat = predict(regression.forest.rep)$predictions, num.trees = num.trees
  )
  causal.forest.weight <- causal_forest(X[rep(train, 2), ], c(Y.a[train], Y.b[train]), c(W.a[train], W.b[train]),
    W.hat = predict(regression.forest.weight)$predictions, sample.weights = sample.weights, num.trees = num.trees
  )
  causal.forest.biased <- causal_forest(X[rep(train, 2), ], c(Y.a[train], Y.b[train]), c(W.a[train], W.b[train]),
    W.hat = predict(regression.forest.biased)$predictions, num.trees = num.trees
  )
  expect_true(mean(z_scores(
    predict(causal.forest.rep, X[test, ], estimate.variance = TRUE),
    predict(causal.forest.weight, X[test, ], estimate.variance = TRUE)
  ) <= 1) >= .5)
  expect_true(mean(predict(causal.forest.rep, X[test, ]) > 100 + predict(causal.forest.biased, X[test, ])) >= .5)
})

test_that("A non-pruned honest causal forest contains trees with empty leafs,
          and a pruned honest causal forest does not contain trees with empty leafs", {
  n <- 100
  p <- 4
  num.trees <- 100
  trees <- 1:num.trees
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- 2 * X[, 1] + X[, 2]
  Y <- W * tau + rnorm(n)

  cf.unpruned <- causal_forest(X, Y, W, honesty = TRUE, honesty.fraction = 0.9,
                               honesty.prune.leaves = FALSE, num.trees = num.trees)
  cf.pruned <- causal_forest(X, Y, W, honesty = TRUE, honesty.fraction = 0.9,
                             honesty.prune.leaves = TRUE, num.trees = num.trees)

  contains_empty_leafs <- function(forest, trees) {
    empty <- lapply(trees, function(t) {
      tree <- get_tree(forest, t)
      printed.tree <- capture.output(print(tree))
      has.empty.leafs <- any(grepl("num_samples: 0", printed.tree))
      has.empty.leafs
    })
    empty
  }

  empty.unpruned <- contains_empty_leafs(cf.unpruned, trees)
  any.unpruned.empty <- any(as.logical(empty.unpruned))

  empty.pruned <- contains_empty_leafs(cf.pruned, trees)
  any.pruned.empty <- any(as.logical(empty.pruned))

  expect_true(any.unpruned.empty)
  expect_true(!any.pruned.empty)
})
