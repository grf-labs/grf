library(grf)

test_that("a simple survival forest workflow works", {
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  failure.time <- -log(runif(n)) * exp(0.1 * X[, 1])
  censor.time <- rexp(n)
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)

  sf <- survival_forest(X, Y, D, num.trees = 50)
  n.failures <- length(sf[["failure.times"]])

  survival.oob <- predict(sf)$predictions
  survival <- predict(sf, X)$predictions
  survival.na <- predict(sf, X, prediction.type = "Nelson-Aalen")$predictions

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

  # Nelson-Aalen estimates of the survival curve are above zero
  expect_true(all(survival.na > 0))

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

  # Predicting at the same custom grid as the training grid gives the same result
  failure.times <- sf$failure.times
  survival.oob.grid <- predict(sf, failure.times = failure.times)$predictions
  survival.grid <- predict(sf, X, failure.times = failure.times)$predictions

  expect_equal(survival.oob.grid, survival.oob)
  expect_equal(survival.grid, survival)
})

test_that("survival forest grid indexing works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  failure.time <- round(exp(0.5 * X[, 1]) * rexp(n), 2)
  censor.time <- round(2 * rexp(n), 2)
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)
  sf <- survival_forest(X, Y, D, num.trees = 250)
  
  # Pick out S(Yi, X) from the estimated survival curve in two identical ways
  sf.pred1 <- predict(sf)
  grid.index1 <- sf$Y.relabeled
  S.hat1 <- sf.pred1$predictions
  S.Y.hat1 <- rep(1, n)
  S.Y.hat1[grid.index1 != 0] <- S.hat1[cbind(1:n, grid.index1)]
  
  sf.pred2 <- predict(sf, failure.times = sort(unique(Y)))
  grid.index2 <- match(Y, sf.pred2$failure.times)
  S.hat2 <- sf.pred2$predictions
  S.Y.hat2 <- S.hat2[cbind(1:n, grid.index2)]
  
  expect_equal(S.Y.hat1, S.Y.hat2)
})

test_that("sample weighted survival prediction is invariant to weight rescaling", {
  # Estimates of the survival function adjusting each sample count
  # by its rescaled sample weight should leave predictions unchanged.
  # (We can not do a check on one forest with some samples duplicated as the splits might be different)
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- 10 * round(1 + rexp(n), 1)
  D <- rbinom(n, 1, 0.5)

  sample.weights <- runif(n)

  sf1 <- survival_forest(X, Y, D,
                         sample.weights = sample.weights,
                         num.trees = 50,
                         seed = 1)

  sf2 <- survival_forest(X, Y, D,
                        sample.weights = runif(1) * sample.weights,
                        num.trees = 50,
                        seed = 1)
  type <- "Kaplan-Meier"
  expect_equal(predict(sf1, prediction.type = type)$predictions,
               predict(sf2, prediction.type = type)$predictions, tolerance = 1e-8)
  expect_equal(predict(sf1, X, prediction.type = type)$predictions,
               predict(sf2, X, prediction.type = type)$predictions, tolerance = 1e-8)
  type <- "Nelson-Aalen"
  expect_equal(predict(sf1, prediction.type = type)$predictions,
               predict(sf2, prediction.type = type)$predictions, tolerance = 1e-8)
  expect_equal(predict(sf1, X, prediction.type = type)$predictions,
               predict(sf2, X, prediction.type = type)$predictions, tolerance = 1e-8)
})

test_that("survival_forest works as expected with missing values", {
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  failure.time <- -log(runif(n)) * exp(0.1 * X[, 1])
  censor.time <- rexp(n)
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)
  nmissing <- 500
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  sf <- survival_forest(X, Y, D, num.trees = 500)

  # MIA with data duplication
  Xl <- X
  Xr <- X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)
  sf.mia <- survival_forest(X.mia, Y, D, num.trees = 500)

  mse.oob.diff <- mean(rowMeans((predict(sf)$predictions - predict(sf.mia)$predictions)^2))
  mse.diff <- mean(rowMeans((predict(sf, X)$predictions - predict(sf.mia, X.mia)$predictions)^2))

  expect_equal(mse.oob.diff, 0, tolerance = 0.001)
  expect_equal(mse.diff, 0, tolerance = 0.001)
})

test_that("survival forest with complete data is ~equal to regression forest", {
  n <- 500
  p <- 5
  X <- matrix(runif(n * p), n, p)
  Y.max <- 1
  failure.time <- pmin(rexp(n) * X[, 1], Y.max)
  censor.time <- 1e6 * runif(n)
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)

  sf <- survival_forest(X, Y, D, num.trees = 500)
  pp.sf <- predict(sf)
  Y.hat.sf <- expected_survival(pp.sf$predictions, pp.sf$failure.times) #integral of the survival function

  rf <- regression_forest(X, Y, num.trees = 500)
  pp.rf <- predict(rf)$predictions

  expect_equal(Y.hat.sf, pp.rf, tolerance = 0.075)
})
