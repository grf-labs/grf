library(grf)

test_that("a simple survival forest workflow works", {
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  latent.time <- -log(runif(n)) * exp(0.1 * X[, 1])
  censor.time <- rexp(n)
  Y <- pmin(latent.time, censor.time)
  delta <- as.integer(latent.time <= censor.time)

  sf <- survival_forest(X, Y, delta, num.trees = 50)
  n.failures <- length(sf[["failure.times"]])

  survival.oob <- predict(sf)$predictions
  survival <- predict(sf, X)$predictions

  # Predictions are monotonically decreasing
  if (n.failures > 1) {
    expect_true(all(survival.oob[, 1:(n.failures - 1)] - survival.oob[, 2:n.failures] >= 0))
    expect_true(all(survival[, 1:(n.failures - 1)] - survival[, 2:n.failures] >= 0))
  }

  # Prediction dimensions are consistent
  expect_equal(c(n, n.failures),
               dim(survival.oob))
  expect_equal(c(n, n.failures),
             dim(survival))

  # A tree with no failures does not split
  sf <- survival_forest(X, Y, rep(0, n), num.trees = 50)
  n.failures <- length(sf[["failure.times"]])

  survival.oob <- predict(sf)$predictions
  survival <- predict(sf, X)$predictions

  # Prediction dimensions are consistent
  expect_equal(c(n, n.failures),
              dim(survival.oob))
  expect_equal(c(n, n.failures),
              dim(survival))

  tree <- get_tree(sf, 1)
  expect_equal(length(tree[["nodes"]]), 1)
})

test_that("sample weighted survival prediction is invariant to weight rescaling", {
  # With Kaplan-Meier estimates of the survival function adjusting each sample count
  # by its rescaled sample weight should leave predictions unchanged.

  # We can not do a check on one forest with some samples duplicated as the splits might be different
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- 10 * round(1 + rexp(n), 1)
  delta <- rbinom(n, 1, 0.5)

  sample.weights <- runif(n)

  sf1 <- survival_forest(X, Y, delta,
                         sample.weights = sample.weights,
                         num.trees = 50,
                         seed = 1)

  sf2 <- survival_forest(X, Y, delta,
                        sample.weights = runif(1) * sample.weights,
                        num.trees = 50,
                        seed = 1)

  expect_true(sum(abs(predict(sf1, X)$pred - predict(sf2, X)$pred)) < 1e-8)
})

test_that("survival_forest works as expected with missing values", {
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  latent.time <- -log(runif(n)) * exp(0.1 * X[, 1])
  censor.time <- rexp(n)
  Y <- pmin(latent.time, censor.time)
  delta <- as.integer(latent.time <= censor.time)
  nmissing <- 500
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  sf <- survival_forest(X, Y, delta, num.trees = 500)

  # MIA with data duplication
  Xl <- X
  Xr <- X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)
  sf.mia <- survival_forest(X.mia, Y, delta, num.trees = 500)

  mse.oob.diff <- mean(rowMeans((predict(sf)$predictions - predict(sf.mia)$predictions)^2))
  mse.diff <- mean(rowMeans((predict(sf, X)$predictions - predict(sf.mia, X.mia)$predictions)^2))

  expect_equal(mse.oob.diff, 0, tol = 0.001)
  expect_equal(mse.diff, 0, tol = 0.001)
})
