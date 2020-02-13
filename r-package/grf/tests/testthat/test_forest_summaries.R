library(grf)

set.seed(1234)

test_that("causal forest calibration is reasonable", {
  n <- 800
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 1], 0) * (W - 0.75) + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.25 + 0.5 * (X[, 1] > 0),
    Y.hat = pmax(X[, 1], 0) * (0.5 * (X[, 1] > 0) - 0.5),
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 1] - 1) <= 0.4)
  expect_true(abs(tc[2, 1] - 1) <= 0.4)
})

test_that("causal forest calibration is reasonable with no average effect", {
  n <- 800
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rnorm(n, 1, 0.5)
  Y <- sign(X[, 1]) * (W - 0.5) + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.5,
    Y.hat = 0,
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 3]) <= 3)
  expect_true(abs(tc[2, 1] - 1) <= 0.3)
})

test_that("causal forest calibration is reasonable with no heterogeneous effect", {
  n <- 800
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 2], 0) + W + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.25 + 0.5 * (X[, 1] > 0),
    Y.hat = 0.25 + 0.5 * (X[, 1] > 0) + pmax(X[, 2], 0),
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 1] - 1) <= 0.3)
  expect_true(abs(tc[2, 3]) <= 4)
})

test_that("causal forest calibration is reasonable with no heterogeneous effect with sample weights and clusters", {
  p <- 4
  K <- 100
  cluster.sizes <- pmax(20, round(40 + 3 * rt(K, df = 3)))
  n <- sum(cluster.sizes)
  clust <- rep(1:K, cluster.sizes)

  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  big <- as.numeric(cluster.sizes[clust] >= median(cluster.sizes))
  Y <- pmax(X[, 2], 0) * big + W + rnorm(n)

  e.cc <- 1 / (1 + exp(-2 * X[, 1]))
  cc <- as.logical(rbinom(n, 1, e.cc))
  ipcc.sample.weights <- 1 / e.cc
  clust.sample.weights <- 1 / cluster.sizes[clust]
  sample.weights <- ipcc.sample.weights * clust.sample.weights

  cf <- causal_forest(X[cc, ], Y[cc], W[cc],
    W.hat = 0.25 + 0.5 * (X[cc, 1] > 0),
    Y.hat = 0.25 + 0.5 * (X[cc, 1] > 0) + pmax(X[cc, 2], 0) * big[cc],
    sample.weights = sample.weights[cc],
    clusters = clust[cc],
    equalize.cluster.weights = FALSE,
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 1] - 1) <= 0.3)
  expect_true(abs(tc[2, 3]) <= 4)
})

test_that("regression forest calibration is reasonable", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- 5 + 5 * sign(X[, 1]) + rnorm(n)

  rf <- regression_forest(X, Y)
  tc <- test_calibration(rf)

  expect_true(abs(tc[1, 1] - 1) <= 0.1)
  expect_true(abs(tc[2, 1] - 1) <= 0.2)
})

test_that("regression forest calibration is reasonable with no heterogeneous effect", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- 5 + rnorm(n)

  rf <- regression_forest(X, Y)
  tc <- test_calibration(rf)

  expect_true(abs(tc[1, 1] - 1) <= 0.1)
  expect_true(abs(tc[2, 3]) <= 4)
})

test_that("causal forest calibration works with clusters", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 2], 0) + W + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.25 + 0.5 * (X[, 1] > 0),
    Y.hat = 0.25 + 0.5 * (X[, 1] > 0) + pmax(X[, 2], 0),
    num.trees = 100
  )
  tc <- test_calibration(cf)

  cf.clust <- cf
  cf.clust$W.orig <- c(cf$W.orig[1:(n / 2)], rep(cf$W.orig[n / 2 + 1:(n / 2)], 10))
  cf.clust$Y.orig <- c(cf$Y.orig[1:(n / 2)], rep(cf$Y.orig[n / 2 + 1:(n / 2)], 10))
  cf.clust$W.hat <- c(cf$W.hat[1:(n / 2)], rep(cf$W.hat[n / 2 + 1:(n / 2)], 10))
  cf.clust$Y.hat <- c(cf$Y.hat[1:(n / 2)], rep(cf$Y.hat[n / 2 + 1:(n / 2)], 10))
  cf.clust$predictions <- c(cf$predictions[1:(n / 2)], rep(cf$predictions[n / 2 + 1:(n / 2)], 10))
  cf.clust$debiased.error <- c(cf$debiased.error[1:(n / 2)], rep(cf$debiased.error[n / 2 + 1:(n / 2)], 10))
  cf.clust$excess.error <- c(cf$excess.error[1:(n / 2)], rep(cf$excess.error[n / 2 + 1:(n / 2)], 10))
  cf.clust$clusters <- c(1:(n / 2), rep(n / 2 + 1:(n / 2), 10))
  cf.clust$equalize.cluster.weights <- TRUE
  tc.clust <- test_calibration(cf.clust)

  expect_equal(tc[, 1:3], tc.clust[, 1:3])
})

test_that("best linear projection is reasonable", {
  n <- 2000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  S <- rbinom(n, 1, 0.5)
  Y <- S * pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  seed = runif(1, 0, .Machine$integer.max)
  forest <- causal_forest(X, Y, W, num.trees = 300, seed = seed)
  forest.shifted.W <- causal_forest(X, Y, W + 1, num.trees = 300, seed = seed)
  forest.2W <- causal_forest(X, Y, W * 2, num.trees = 300, seed = seed)

  # hard-code true regression coefficients
  beta0.true <- 0.2 # approximation to mean(pmax(Z, 0))/2
  beta1.true <- 0.25 # cov(pmax(Z, 0), Z)/2

  blp.all <- best_linear_projection(forest, X[,1:2])
  expect_equal(blp.all[1,1], beta0.true, tol = 0.1)
  expect_equal(blp.all[2,1], beta1.true, tol = 0.1)
  expect_equal(blp.all[3,1], 0, tol = 0.1)

  blp.subset <- best_linear_projection(forest, X[,1:2], subset = (S == 1))
  expect_equal(blp.subset[1,1], beta0.true * 2, tol = 0.2)
  expect_equal(blp.subset[2,1], beta1.true * 2, tol = 0.2)
  expect_equal(blp.subset[3,1], 0, tol = 0.2)

  std.errs <- c((blp.all[1,1] - beta0.true) / blp.all[1,2],
               (blp.all[2,1] - beta1.true) / blp.all[2,2],
               (blp.all[3,1] - 0) / blp.all[3,2],
               (blp.subset[1,1] - beta0.true * 2) / blp.subset[1,2],
               (blp.subset[2,1] - beta1.true * 2) / blp.subset[2,2],
               (blp.subset[3,1] - 0) / blp.subset[3,2])

  expect_lt(mean(abs(std.errs)), 2.5)
  expect_gt(mean(abs(std.errs)), 0.25)

  blp.shifted <- best_linear_projection(forest.shifted.W, X[,1:2])
  expect_equal(blp.all[, "Estimate"], blp.shifted[, "Estimate"], tol = 0.05)
  expect_equal(blp.all[, "Std. Error"], blp.shifted[, "Std. Error"], tol = 0.05)

  blp.2W <- best_linear_projection(forest.2W, X[,1:2])
  expect_equal(blp.all[, "Estimate"]/2, blp.2W[, "Estimate"], tol = 0.05)
})

test_that("best linear projection works with edge case input types", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  forest <- causal_forest(X, Y, W, num.trees = 50)
  forest.clustered <- causal_forest(X, Y, W, num.trees = 50, clusters = c(rep(1, 100), rep(2, 100)))

  # a single covariate
  blp.single.covariate <- best_linear_projection(forest, X[, 1])
  # a forest trained with clusters, and subset equal to only one of the clusters
  blp.single.cluster <- best_linear_projection(forest, X[, 1], subset = which(forest.clustered$clusters == 1))
  expect_equal(1, 1)
})
