library(grf)

seed <- 1000
set.seed(1000)

test_that("Merged regression forest attributes are sensible", {
  # Train regression forests
  seed <- 1000
  set.seed(seed)
  n = 50; p = 2
  X = matrix(rnorm(n*p), n, p)
  Y = X[,1] * rnorm(n)
  W = X[,2] > 0
  r.forest1 = regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10, seed = seed)
  r.forest2 = regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10, seed = seed)

  # Join the forests together.
  big.rf = merge_forests(list(r.forest1, r.forest2))

  # Result is also a regression_forest of the same class
  expect_true(is(big.rf, "grf"))
  expect_equal(r.forest1[["_num_trees"]] + r.forest2[["_num_trees"]], big.rf[["_num_trees"]])
  expect_equal(class(r.forest1), class(big.rf))
})

test_that("Merged causal forest attributes are sensible", {
  # Train causal forests
  seed <- 1000
  set.seed(seed)
  n = 150; p = 3
  X = matrix(rnorm(n*p), n, p)
  Y = X[,1] * rnorm(n)
  W = X[,2] > 0
  c.forest1 = causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 25, seed = seed)
  c.forest2 = causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 25, seed = seed)

  # Join the forests together.
  big.rf = merge_forests(list(c.forest1, c.forest2))

  # Result is also a causal forest of the same class
  expect_true(big.rf[["_num_trees"]] == (c.forest1[["_num_trees"]] + c.forest2[["_num_trees"]]))
  expect_true(is(big.rf, "grf"))
  expect_equal(class(c.forest1), class(big.rf))

  expect_equal(c.forest1$Y.hat, big.rf$Y.hat)
  expect_equal(c.forest1$W.hat, big.rf$W.hat)
})


test_that("Merged causal forests give reasonable predictions", {
  seed <- 1000
  set.seed(seed)
  n = 200; p = 15
  X = matrix(rnorm(n*p), n, p)
  Y = X[,1] * X[,3] * X[,4] * rnorm(n) + X[,6] * X[,7] * rnorm(n)
  W = X[,2] * X[,5] > 0

  # Train a causal forest.
  num.trees <- 100
  num.forests <- 20
  c.forests <- lapply(seq(num.forests), function(i) causal_forest(X, Y, W, num.trees = num.trees, seed = i))
  big.rf = merge_forests(c.forests)

  c.forest1 <- c.forests[[1]]
  preds = predict(c.forest1)
  excess.error = mean(preds$excess.error)
  error = mean(preds$debiased.error) + excess.error

  big.preds = predict(big.rf)
  big.excess.error = mean(big.preds$excess.error)
  big.error = mean(big.preds$debiased.error) + big.excess.error

  expect_lt(big.error, error)
  expect_lt(big.excess.error, excess.error)
})

test_that("Incompatible forests are not mergeable", {
  # Train causal forests
  seed <- 1000
  set.seed(seed)
  n = 50; p = 2
  X = matrix(rnorm(n*p), n, p)
  Y = X[,1] * rnorm(n)
  W = X[,2] > 0

  c.forest1 = causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 10,
                            ci.group.size = 3, seed = seed)
  c.forest2 = causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 10,
                            ci.group.size = 4, seed = seed)
  r.forest1 = regression_forest(X, Y, ci.group.size = 3, seed = seed)

  #  Empty input
  expect_error(merge_forests(list()))
  # Incompatible classes
  expect_error(merge_forests(list(c.forest1, r.forest1)))
  # Incompatible ci.group.size
  expect_error(merge_forests(list(c.forest1, c.forest2)))
})
