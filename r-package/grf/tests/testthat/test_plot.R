test_that("basic quantile forest plotting is successful", {
  p <- 4
  n <- 50
  i <- 2
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))
  D <- data.frame(X = X, Y = Y)
  q.forest <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), num.trees = 50)
  q.tree <- get_tree(q.forest, 1)
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    capture_output(plot(q.tree))
  }
  expect_true(TRUE)
})

test_that("basic regression forest plotting is successful", {
  p <- 4
  n <- 50
  i <- 2
  X <- matrix(2 * runif(n * p) - 1, n, p)
  colnames(X) <- c("age", "income", "weight", "height")
  Y <- rnorm(n) * (1 + (X[, i] > 0))
  r.forest <- regression_forest(X, Y)
  r.tree <- get_tree(r.forest, 1)
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    capture_output(plot(r.tree))
  }
  expect_true(TRUE)
})

test_that("large regression forest plotting is successful", {
  p <- 4
  n <- 500
  i <- 2
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))
  r.forest <- regression_forest(X, Y)
  r.tree <- get_tree(r.forest, 1)
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    capture_output(plot(r.tree))
  }
  expect_true(TRUE)
})
