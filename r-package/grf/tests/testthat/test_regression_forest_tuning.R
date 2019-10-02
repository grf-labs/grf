test_that("regression forest tuning decreases prediction error", {
  n <- 5000
  p <- 2

  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- (X[, 1] > 0) + rnorm(n)
  X.test <- matrix(2 * runif(n * p) - 1, n, p)
  truth <- (X.test[, 1] > 0)

  forest <- regression_forest(X, Y, num.trees = 400, tune.parameters = "none")
  preds <- predict(forest, X.test)
  error <- mean((preds$predictions - truth)^2)

  tuned.forest <- regression_forest(X, Y, num.trees = 400, tune.parameters = "all")
  tuned.preds <- predict(tuned.forest, X.test)
  tuned.error <- mean((tuned.preds$predictions - truth)^2)

  expect_true(tuned.error < error * 0.75)
})

test_that("regression forest tuning tunes each parameter", {
  n <- 100
  p <- 2

  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- (X[, 1] > 0) + rnorm(n)
  tunable.params <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
                    "honesty.prune.leaves", "alpha", "imbalance.penalty")
  for (param in tunable.params) {
    tuned.forest <- regression_forest(X, Y, num.trees = 100, tune.parameters = param)
    expect_true(param == names(tuned.forest$tuning.output$params))
  }
})
