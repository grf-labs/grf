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

  expect_lt(tuned.error, error * 0.75)
})

test_that("tuning with a seed ensures reproducible results without disturbing global seed", {
  n <- 200
  p <- 2

  X <- matrix(runif(n * p), n, p)
  Y <- (X[, 1] > 0) + rnorm(n)

  global.seed.before <- .Random.seed
  tuned.forest.1 <- regression_forest(X, Y, num.trees = 150, tune.parameters = "all", seed = 1)
  tuned.preds.1 <- predict(tuned.forest.1)
  tuned.forest.2 <- regression_forest(X, Y, num.trees = 150, tune.parameters = "all", seed = 1)
  tuned.preds.2 <- predict(tuned.forest.2)
  global.seed.after <- .Random.seed

  expect_equal(tuned.preds.1, tuned.preds.2)
  expect_equal(global.seed.before, global.seed.after)
})
