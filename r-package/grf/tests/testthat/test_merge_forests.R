library(grf)

test_that("Merged regression forest attributes are sensible", {
  # Train regression forests
  n <- 50
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)
  W <- X[, 2] > 0
  r.forest1 <- regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10)
  r.forest2 <- regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10)

  # Join the forests together.
  big.rf <- merge_forests(list(r.forest1, r.forest2))

  # Result is also a regression_forest of the same class
  expect_true(methods::is(big.rf, "grf"))
  expect_equal(r.forest1[["_num_trees"]] + r.forest2[["_num_trees"]], big.rf[["_num_trees"]])
  expect_equal(class(r.forest1), class(big.rf))
})

test_that("Merged causal forest attributes are sensible", {
  # Train causal forests
  n <- 50
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)
  W <- X[, 2] > 0
  c.forest1 <- causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 10)
  c.forest2 <- causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 10)

  # Join the forests together.
  big.rf <- merge_forests(list(c.forest1, c.forest2))

  # Result is also a causal forest of the same class
  expect_true(big.rf[["_num_trees"]] == (c.forest1[["_num_trees"]] + c.forest2[["_num_trees"]]))
  expect_true(methods::is(big.rf, "grf"))
  expect_equal(class(c.forest1), class(big.rf))

  expect_equal(c.forest1$Y.hat, big.rf$Y.hat)
  expect_equal(c.forest1$W.hat, big.rf$W.hat)
})


test_that("Merged causal forests give reasonable predictions", {
  n <- 1000
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  W <- X[, 2] > 0
  Tau <- X[, 3]
  Y <- X[, 1] * rnorm(n) + W * X[, 2] + rnorm(n)
  num.trees <- 25

  # Train a causal forest.
  c.forest1 <- causal_forest(X, Y, W, num.trees = num.trees, alpha = 0, min.node.size = 1)

  # Train another sequence of forests of equal size, then merge them.
  c.forests <- lapply(seq(100), function(x) causal_forest(X, Y, W, alpha = 0, num.trees = num.trees, min.node.size = 1))
  big.rf <- merge_forests(c.forests)

  # The merged forest should have smaller error
  preds <- predict(c.forest1)
  excess.error <- mean(preds$excess.error, na.rm = TRUE)
  error <- mean((Tau - preds$predictions)^2, na.rm = TRUE)

  big.preds <- predict(big.rf)
  big.excess.error <- mean(big.preds$excess.error)
  big.error <- mean((Tau - big.preds$predictions)^2, na.rm = TRUE)

  expect_lt(big.error / error, 0.95)
  expect_lt(big.excess.error / excess.error, 0.5)
})

test_that("Incompatible forests are not mergeable", {
  # Train causal forests
  n <- 50
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)
  W <- X[, 2] > 0

  c.forest1 <- causal_forest(X, Y, W,
    compute.oob.predictions = FALSE, num.trees = 10,
    ci.group.size = 3
  )
  c.forest2 <- causal_forest(X, Y, W,
    compute.oob.predictions = FALSE, num.trees = 10,
    ci.group.size = 4
  )
  r.forest1 <- regression_forest(X, Y, ci.group.size = 3)

  #  Empty input
  expect_error(merge_forests(list()))
  # Incompatible classes
  expect_error(merge_forests(list(c.forest1, r.forest1)))
  # Incompatible ci.group.size
  expect_error(merge_forests(list(c.forest1, c.forest2)))

  # Forests trained on data sets of different sizes are not mergeable
  r.forest1 <- regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10)
  r.forest2 <- regression_forest(X[, 1, drop = FALSE], Y, compute.oob.predictions = FALSE, num.trees = 10)
  expect_error(merge_forests(list(r.forest1, r.forest2)))

  r.forest1 <- regression_forest(X[1:20, ], Y[1:20], compute.oob.predictions = FALSE, num.trees = 10)
  r.forest1 <- regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10)
  expect_error(merge_forests(list(r.forest1, r.forest2)))
})
