test_that("instrumental forest tuning decreases prediction error", {
  n <- 1500; p <- 5
  X <- matrix(rbinom(n * p, 1, 0.5), n, p)
  Z <- rbinom(n, 1, 0.5)
  Q <- rbinom(n, 1, 0.5)
  T <- Q * Z
  eps <- rnorm(n)
  TAU <- X[, 1] / 2
  Y <- rowSums(X[, 1:3]) + TAU * T + Q + eps
  Y.hat <- predict(regression_forest(X, Y, num.trees = 500))$predictions
  W.hat <- predict(regression_forest(X, T, num.trees = 500))$predictions
  Z.hat <- predict(regression_forest(X, Z, num.trees = 500))$predictions

  iv.forest <- instrumental_forest(X, Y, T, Z, Y.hat, W.hat, Z.hat, tune.parameters = "none")
  error <- mean((TAU - predict(iv.forest)$predictions)^2)

  tuned.iv.forest <- instrumental_forest(X, Y, T, Z, Y.hat, W.hat, Z.hat, tune.parameters = "all")
  tuned.error <- mean((TAU - predict(tuned.iv.forest)$predictions)^2)

  expect_lt(tuned.error, error * 0.8)
})
