library(grf)

test_that("custom forests behave as expected", {
  p <- 40
  n <- 500

  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, 1] > 0))
  D <- data.frame(X = X, Y = Y)

  forest <- custom_forest(X, Y)
  predictions <- predict(forest, X)
  expect_equal(0.0, sum(predictions))
})

test_that("custom forest prediction is of correct size", {
  p <- 40
  n <- 500
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, 1] > 0))

  n_test <- 50
  X_test <- matrix(2 * runif(n_test * p) - 1, n_test, p)

  forest <- custom_forest(X, Y)
  predictions <- predict(forest, X_test)
  expect_length(predictions, n_test)
})
