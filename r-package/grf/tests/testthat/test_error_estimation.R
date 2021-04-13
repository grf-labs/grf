library(grf)

set.seed(1234)

test_that("regression error estimates are reasonable", {
  p <- 3
  n <- 2000
  sigma <- 1

  X <- matrix(2 * runif(n * p) - 1, n, p)
  MU <- 0.1 * (X[, 1] > 0)
  Y <- MU + sigma * rnorm(n)

  output.5.reps <- replicate(20, {
    forest.5 <- regression_forest(X, Y, num.trees = 5, ci.group.size = 1, min.node.size = 1)
    pred.5 <- predict(forest.5)
    c(
      mean(pred.5$debiased.error, na.rm = TRUE),
      mean((pred.5$predictions - Y)^2, na.rm = TRUE)
    )
  })
  err.debiased.5 <- mean(output.5.reps[1, ])
  mse.5 <- mean(output.5.reps[2, ])

  forest.200 <- regression_forest(X, Y, num.trees = 200, ci.group.size = 1)
  pred.200 <- predict(forest.200)
  err.debiased.200 <- mean(pred.200$debiased.error, na.rm = TRUE)
  mse.200 <- mean((pred.200$predictions - Y)^2, na.rm = TRUE)

  expect_equal(err.debiased.200, mse.200, tolerance = 0.01 * sigma^2)
  expect_equal(err.debiased.5, err.debiased.200, tolerance = 0.025 * sigma^2)

  expect_gte(mse.5 - mse.200, sigma^2 / 10)
})

test_that("causal error estimates are reasonable", {
  p <- 3
  n <- 2000
  sigma <- 0.1

  X <- matrix(2 * runif(n * p) - 1, n, p)
  W <- rbinom(n, 1, 0.1)
  TAU <- (X[, 1] > 0)
  Y <- 2 * TAU * (W - 1 / 2) + sigma * rnorm(n)

  W.forest <- regression_forest(X, W, num.trees = 500, sample.fraction = 0.2)
  W.hat <- predict(W.forest)$predictions

  Y.forest <- regression_forest(X, Y, num.trees = 500, sample.fraction = 0.2)
  Y.hat <- predict(Y.forest)$predictions

  Y.resid <- Y - Y.hat
  W.resid <- W - W.hat

  output.10.reps <- replicate(20, {
    cf.10 <- causal_forest(X, Y, W,
      Y.hat = Y.hat, W.hat = W.hat,
      num.trees = 10, min.node.size = 1, stabilize.splits = TRUE
    )
    tau.hat.10 <- predict(cf.10)
    c(
      mean((Y.resid - tau.hat.10$predictions * W.resid)^2, na.rm = TRUE),
      mean(tau.hat.10$debiased.error, na.rm = TRUE)
    )
  })
  raw.10 <- mean(output.10.reps[1, ])
  err.10 <- mean(output.10.reps[2, ])

  output.20.reps <- replicate(10, {
    cf.20 <- causal_forest(X, Y, W,
      Y.hat = Y.hat, W.hat = W.hat,
      num.trees = 20, min.node.size = 1, stabilize.splits = TRUE
    )
    tau.hat.20 <- predict(cf.20)
    c(
      mean((Y.resid - tau.hat.20$predictions * W.resid)^2, na.rm = TRUE),
      mean(tau.hat.20$debiased.error, na.rm = TRUE)
    )
  })
  raw.20 <- mean(output.20.reps[1, ])
  err.20 <- mean(output.20.reps[2, ])

  cf.400 <- causal_forest(X, Y, W,
    Y.hat = Y.hat, W.hat = W.hat,
    num.trees = 400, min.node.size = 1, stabilize.splits = TRUE
  )
  tau.hat.400 <- predict(cf.400)
  raw.400 <- mean((Y.resid - tau.hat.400$predictions * W.resid)^2)
  err.400 <- mean(tau.hat.400$debiased.error, na.rm = TRUE)

  expect_equal(err.400, raw.400, tolerance = 0.1 * sigma^2)
  expect_equal(err.10, err.400, tolerance = 1.5 * sigma^2)
  expect_equal(err.20, err.400, tolerance = 1.0 * sigma^2)
  expect_gt(raw.10 - err.400, sigma^2)
  expect_lt(err.10 - err.400, sigma^2)
})
