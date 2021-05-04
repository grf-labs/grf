test_that("basic printing is successful", {
  p <- 4
  n <- 50
  i <- 2
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))
  D <- data.frame(X = X, Y = Y)
  q.forest <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), num.trees = 50)
  capture_output(print(q.forest))
  expect_true(TRUE)
})

test_that("printing a forest with one regressor is successful", {
  n <- 50
  p <- 1
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)
  r.forest <- regression_forest(X, Y)
  capture_output(print(r.forest))
  expect_true(TRUE)
})

test_that("basic tree printing is successful", {
  p <- 4
  n <- 50
  i <- 2
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))
  D <- data.frame(X = X, Y = Y)
  q.forest <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), num.trees = 50)
  q.tree <- get_tree(q.forest, 1)
  capture_output(print(q.tree))
  expect_true(TRUE)
})

test_that("tuning output printing is successful", {
  p <- 4
  n <- 100
  X <- matrix(runif(n * p), n, p)
  Y <- runif(n)
  W <- rbinom(n, 1, 0.5)
  tuned.forest <- causal_forest(X, Y, W, num.trees = 100, tune.parameters = "all")
  capture_output(print(tuned.forest$tuning.output))
  expect_true(TRUE)
})
