library(grf)

set.seed(3141)

test_that("causal forests can split on the last parameter", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * (X[, 1] + X[, 6]) + rnorm(n)

  forest <- causal_forest(X, Y, W, compute.oob.predictions = FALSE)
  split.freq <- split_frequencies(forest, 10)

  expect_gt(sum(split.freq[, 6]), 0)
})

test_that("causal forests have reasonable split frequencies", {
  n <- 100
  p <- 7
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.2)
  Y <- 1000 * (X[, p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase imbalance.penalty to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc <- causal_forest(X, Y, W, mtry = p, imbalance.penalty = 1.0, stabilize.splits = TRUE, min.node.size = 2)
  split.freq <- split_frequencies(ccc, 4)
  expect_gt(split.freq[1, p] / sum(split.freq[1, ]), 2 / 3)
})

test_that("causal forests without stable splitting have reasonable split frequencies", {
  n <- 100
  p <- 7
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.2)
  Y <- 1000 * (X[, p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase imbalance.penalty to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc <- causal_forest(X, Y, W, mtry = p, imbalance.penalty = 1.0, stabilize.splits = FALSE, min.node.size = 2)
  split.freq <- split_frequencies(ccc, 4)
  expect_gt(split.freq[1, p] / sum(split.freq[1, ]), 2 / 3)
})

test_that("causal forests with a positive imbalance.penalty have reasonable tree depths", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(p = 0.5, size = 1, n = n)
  Y <- 0.5 * X[, 1] * (2 * W - 1) + 0.1 * rnorm(n)

  forest <- causal_forest(X, Y, W, imbalance.penalty = 0.001)
  split.freq <- split_frequencies(forest)
  expect_gt(sum(split.freq[3, ]), 0)
})

test_that("causal forests with a very small imbalance.penalty behave similarly to unpenalized forests.", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- X[, 1] * (2 * W - 1) + 0.1 * rnorm(n)

  forest <- causal_forest(X, Y, W, imbalance.penalty = 0.0)
  forest.large.penalty <- causal_forest(X, Y, W, imbalance.penalty = 100.0)
  forest.small.penalty <- causal_forest(X, Y, W, imbalance.penalty = 1e-6)

  diff.large.penalty <- abs(forest.large.penalty$debiased.error - forest$debiased.error)
  diff.small.penalty <- abs(forest.small.penalty$debiased.error - forest$debiased.error)
  expect_lt(mean(diff.small.penalty), 0.10 * mean(diff.large.penalty))
})

test_that("causal forests behave reasonably with a low treatment probability", {
  n <- 1000
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.1)
  tau <- 0.1
  Y <- X[, 1] + X[, 2] + tau * W + rnorm(n)

  forest <- causal_forest(X, Y, W, stabilize.splits = TRUE)
  tau.hat <- predict(forest)$predictions
  expect_lt(sqrt(mean((tau.hat - tau)^2)), 0.20)
})

test_that("causal forests behave reasonably with small sample size", {
  p <- 5
  n <- 50
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- 100 * (X[, 1] > 0)
  Y <- tau * (W - 0.5) + rnorm(n)
  forest <- causal_forest(X, Y, W,
    stabilize.splits = TRUE,
    min.node.size = 1, mtry = p,
    compute.oob.predictions = FALSE
  )

  X.test <- matrix(rnorm(n * p), n, p)
  tau.test <- 100 * (X.test[, 1] > 0)
  tau.hat <- predict(forest, X.test)$predictions
  expect_lt(sqrt(mean((tau.hat - tau.test)^2)) / 100, 1 / 3)
})

test_that("local linear causal forests work in a simple case", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  TAU <- 2 * X[, 1] + X[, 2]
  Y <- W * TAU + rnorm(n)

  forest <- causal_forest(X, Y, W, num.trees = 400)
  preds.ll <- predict(forest, X, linear.correction.variables = 1:2, ll.lambda = 0.01)
  error.ll <- mean((preds.ll$predictions - TAU)^2)

  preds.rf <- predict(forest, X)
  error.rf <- mean((preds.rf$predictions - TAU)^2)

  expect_lt(error.ll, 0.8 * error.rf)
})

test_that("local linear causal forests with large lambda are equivalent to causal forests", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- 2 * X[, 1] + X[, 2]
  Y <- W * tau + rnorm(n)

  forest <- causal_forest(X, Y, W, num.trees = 400)
  preds.ll <- predict(forest, X, linear.correction.variables = 1:ncol(X), ll.lambda = 1e5)$predictions
  preds.cf <- predict(forest)$predictions

  expect_lt(mean((preds.ll - preds.cf)^2), 0.02)
})

test_that("causal forest predictions and variance estimates are invariant to scaling of the sample weights.", {
  n <- 200
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-2 * X[, 1] + 3 * X[, 2]))
  W <- rbinom(n, 1, e)
  tau <- 2 * X[, 1] + X[, 2]
  Y <- W * tau + rnorm(n)
  e.cc <- 1 / (1 + exp(-2 * X[, 1]))
  sample.weights <- 1 / e.cc

  # The multiple is a power of 2 to avoid rounding errors allowing for exact comparison
  # between two forest with the same seed.
  forest.1 <- causal_forest(X, Y, W, sample.weights = sample.weights, seed = 1)
  forest.2 <- causal_forest(X, Y, W, sample.weights = 64 * sample.weights, seed = 1)
  pred.1 <- predict(forest.1, estimate.variance = TRUE)
  pred.2 <- predict(forest.2, estimate.variance = TRUE)

  expect_equal(pred.1$predictions, pred.2$predictions, tolerance = 1e-10)
  expect_equal(pred.1$variance.estimates, pred.2$variance.estimates, tolerance = 1e-10)
  expect_equal(pred.1$debiased.error, pred.2$debiased.error, tolerance = 1e-10)
})

test_that("sample weighted causal forest is estimated with kernel weights `forest.weights * sample.weights`", {
  n <- 500
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-X[, 1] + 2 * X[, 2]))
  W <- rbinom(n, 1, e)
  tau <- 2 * (X[, 1] > 0 & X[, 5] > 0) -
    0.5 * (X[, 2] > 0) - 0.5 * (X[, 3] > 0) - 0.5 * (X[, 4] > 0)
  Y <- W * tau + rnorm(n)
  e.cc <- 1 - 0.9 * (X[, 1] > 0 & X[, 5] > 0)
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc
  cf <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0, sample.weights = sample.weights, num.trees = 250)

  x1 <- X[1, , drop = F]
  theta1 <- predict(cf, x1)$predictions
  alpha1 <- get_forest_weights(cf, x1)[1, ]
  theta1.lm <- lm(Y ~ W, weights = alpha1 * sample.weights)

  expect_equal(theta1, theta1.lm$coefficients[[2]], tolerance = 1e-10)
})

test_that("IPCC weighting in the training of a causal forest with missing data improves its complete-data MSE.", {
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-X[, 1] + 2 * X[, 2]))
  W <- rbinom(n, 1, e)
  tau <- 2 * (X[, 1] > 0 & X[, 5] > 0) -
      0.5 * (X[, 2] > 0) - 0.5 * (X[, 3] > 0) - 0.5 * (X[, 4] > 0)
  Y <- W * tau + rnorm(n)

  e.cc <- 1 - 0.9 * (X[, 1] > 0 & X[, 5] > 0)
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc

  num.trees <- 500
  mse <- function(f) {
      tau.hat <- rep(NA, n)
      tau.hat[cc] <- predict(f)$predictions
      tau.hat[!cc] <- predict(f, X[!cc, ])$predictions
      mean((tau.hat - tau)^2)
  }

  forest <- causal_forest(X[cc, ], Y[cc], W[cc], num.trees = num.trees)
  weighted.forest <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights = sample.weights[cc], num.trees = num.trees)
  expect_lt(mse(weighted.forest) / mse(forest), .9)
})

test_that("sample weighted causal forest gives correct coverage", {
  n <- 2500
  p <- 5
  pA <- 0.2
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, pA)
  gamma <- A / pA + (1 - A) / (1 - pA)
  tau.true <- 1 / (1 + exp(-3 * X[,1]))
  Y <- 2 * tau.true * W * A + rnorm(n)

  cf <- causal_forest(X, Y, W, W.hat = 0.5, sample.weights = gamma)
  preds <- predict(cf, estimate.variance = TRUE)
  zstat <- (preds$predictions - tau.true) / sqrt(preds$variance.estimates)
  coverage <- mean(abs(zstat) < 1.96)
  mse <- mean((preds$predictions - tau.true)^2)

  cf.noweight <- causal_forest(X, Y, W, W.hat = 0.5)
  preds.noweight <- predict(cf.noweight, estimate.variance = TRUE)
  zstat.noweight <- (preds.noweight$predictions - tau.true) / sqrt(preds.noweight$variance.estimates)
  coverage.noweight <- mean(abs(zstat.noweight) < 1.96)
  mse.noweight <- mean((preds.noweight$predictions - tau.true)^2)

  expect_lt(mse / mse.noweight, 0.3)
  expect_gt(coverage, 0.75)
  expect_lt(coverage.noweight, 0.55)
})

test_that("Weighting is roughly equivalent to replication of samples", {
  n <- 500
  p <- 2
  num.trees <- 1000

  X <- matrix(rnorm(n * p), n, p)
  e.a <- 1 / (1 + exp(+1 * X[, 1] + 1 * X[, 2]))
  e.b <- 1 / (1 + exp(-1 * X[, 1] - 1 * X[, 2]))
  tau.a <- +1000 + 2 * X[, 1] + X[, 2]
  tau.b <- -1000 + X[, 1] - X[, 2]
  W.a <- rbinom(n, 1, e.a)
  W.b <- rbinom(n, 1, e.b)
  Y.a <- W.a * tau.a + rnorm(n)
  Y.b <- W.b * tau.b + rnorm(n)

  X.rep <- X[rep(1:(n / 2), 5), ]
  W.rep <- c(W.a[rep(1:(n / 2), 4)], W.b[1:(n / 2)])
  Y.rep <- c(Y.a[rep(1:(n / 2), 4)], Y.b[1:(n / 2)])
  sample.weights <- c(rep(4, n / 2), rep(1, n / 2))
  train <- 1:(n / 2)
  test <- (n / 2 + 1):n

  ## regression forest for propensity score
  regression.forest.rep <- regression_forest(X.rep, W.rep, num.trees = num.trees)
  regression.forest.weight <- regression_forest(X[rep(train, 2), ], c(W.a[train], W.b[train]),
    sample.weights = sample.weights, num.trees = num.trees
  )
  regression.forest.biased <- regression_forest(X[rep(train, 2), ], c(W.a[train], W.b[train]),
    num.trees = num.trees
  )

  z_scores <- function(a, b) {
    # conservative. use a variance upper bound
    abs(a$predictions - b$predictions) / sqrt(2 * (a$variance + b$variance))
  }

  expect_gte(mean(z_scores(predict(regression.forest.rep, X[test, ], estimate.variance = TRUE),
                           predict(regression.forest.weight, X[test, ], estimate.variance = TRUE)) <= 1),
            0.5)
  expect_gte(mean(z_scores(predict(regression.forest.rep, X[test, ], estimate.variance = TRUE),
                           predict(regression.forest.biased, X[test, ], estimate.variance = TRUE)) >= 1),
            0.5)
  ## causal forest
  causal.forest.rep <- causal_forest(X.rep, Y.rep, W.rep,
    W.hat = predict(regression.forest.rep)$predictions, num.trees = num.trees
  )
  causal.forest.weight <- causal_forest(X[rep(train, 2), ], c(Y.a[train], Y.b[train]), c(W.a[train], W.b[train]),
    W.hat = predict(regression.forest.weight)$predictions, sample.weights = sample.weights, num.trees = num.trees
  )
  causal.forest.biased <- causal_forest(X[rep(train, 2), ], c(Y.a[train], Y.b[train]), c(W.a[train], W.b[train]),
    W.hat = predict(regression.forest.biased)$predictions, num.trees = num.trees
  )
  expect_gte(mean(z_scores(predict(causal.forest.rep, X[test, ], estimate.variance = TRUE),
                           predict(causal.forest.weight, X[test, ], estimate.variance = TRUE)) <= 1),
             0.5)
  expect_gte(mean(predict(causal.forest.rep, X[test, ]) > 100 + predict(causal.forest.biased, X[test, ])), 0.5)
})

test_that("A non-pruned honest causal forest contains trees with empty leafs,
          and a pruned honest causal forest does not contain trees with empty leafs", {
  n <- 100
  p <- 4
  num.trees <- 100
  trees <- 1:num.trees
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- 2 * X[, 1] + X[, 2]
  Y <- W * tau + rnorm(n)

  cf.unpruned <- causal_forest(X, Y, W, honesty = TRUE, honesty.fraction = 0.9,
                               honesty.prune.leaves = FALSE, num.trees = num.trees)
  cf.pruned <- causal_forest(X, Y, W, honesty = TRUE, honesty.fraction = 0.9,
                             honesty.prune.leaves = TRUE, num.trees = num.trees)

  contains_empty_leafs <- function(forest, trees) {
    empty <- lapply(trees, function(t) {
      tree <- get_tree(forest, t)
      printed.tree <- capture.output(print(tree))
      has.empty.leafs <- any(grepl("num_samples: 0", printed.tree))
      has.empty.leafs
    })
    empty
  }

  empty.unpruned <- contains_empty_leafs(cf.unpruned, trees)
  any.unpruned.empty <- any(as.logical(empty.unpruned))

  empty.pruned <- contains_empty_leafs(cf.pruned, trees)
  any.pruned.empty <- any(as.logical(empty.pruned))

  expect_true(any.unpruned.empty)
  expect_true(!any.pruned.empty)
})

test_that("causal_forest works as expected with missing values", {
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

  nmissing <- 500
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  # MIA with data duplication
  Xl <- X
  Xr <- X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)

  X.test <- matrix(rnorm(n * p), n, p)
  X.test[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN
  Xlt <- X.test
  Xrt <- X.test
  Xlt[is.nan(Xlt)] <- -1e9
  Xrt[is.nan(Xrt)] <- 1e9
  X.mia.test <- cbind(Xlt, Xrt)

  rf.mia <- causal_forest(X.mia, Y, W, 0, 0, seed = 123)
  rf <- causal_forest(X, Y, W, 0, 0, seed = 123)

  mse.oob.diff <- mean((predict(rf.mia)$predictions - predict(rf)$predictions)^2)
  mse.diff <- mean((predict(rf.mia, X.mia)$predictions - predict(rf, X)$predictions)^2)
  mse.test.diff <- mean((predict(rf.mia, X.mia.test)$predictions - predict(rf, X.test)$predictions)^2)

  diff.mse.oob <- mean((predict(rf.mia)$predictions - Y)^2) - mean((predict(rf)$predictions - Y)^2)
  diff.mse <- mean((predict(rf.mia, X.mia)$predictions - Y)^2) - mean((predict(rf, X)$predictions - Y)^2)
  diff.mse.test <- mean((predict(rf.mia, X.mia.test)$predictions - Y)^2) - mean((predict(rf, X.test)$predictions - Y)^2)

  expect_equal(mse.oob.diff, 0, tolerance = 0.005)
  expect_equal(mse.diff, 0, tolerance = 0.005)
  expect_equal(mse.test.diff, 0, tolerance = 0.005)

  expect_equal(diff.mse.oob, 0, tolerance = 0.05)
  expect_equal(diff.mse, 0, tolerance = 0.05)
  expect_equal(diff.mse.test, 0, tolerance = 0.05)

  # All NaNs
  X[, ] <- NaN

  Xl <- X
  Xr <- X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)

  rf.mia <- causal_forest(X.mia, Y, W, 0, 0, seed = 123)
  rf <- causal_forest(X, Y, W, 0, 0, seed = 123)
  mse.oob.diff.allnan <- mean((predict(rf.mia)$predictions - predict(rf)$predictions)^2)
  expect_equal(mse.oob.diff.allnan, 0, tolerance = 0.0001)
})

test_that("a causal forest workflow with missing values works as expected", {
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

  nmissing <- 2000
  idx.missing <- cbind(sample(1:n, nmissing, replace = TRUE),
                       sample(1:p, nmissing, replace = TRUE))
  X[idx.missing] <- NA

  forest <- causal_forest(X, Y, W, tune.parameters = "all", num.trees = 500)
  tau.hat <- predict(forest)$predictions
  high.effect <- tau.hat > median(tau.hat)

  cal <- test_calibration(forest)
  varimp <- variable_importance(forest)

  ate1 <- average_treatment_effect(forest, subset = high.effect)
  ate2 <- average_treatment_effect(forest, subset = !high.effect)
  ate3 <- average_treatment_effect(forest, subset = is.na(X[, 1]))
  ate4 <- average_treatment_effect(forest, subset = !is.na(X[, 1]))
  ate5 <- average_treatment_effect(forest, subset = complete.cases(X))
  ate6 <- average_treatment_effect(forest, subset = !complete.cases(X))

  blp1 <- best_linear_projection(forest, A = X, subset = complete.cases(X))
  blp2 <- best_linear_projection(forest, A = X[, 1], subset = !is.na(X[, 1]))
  blp3 <- best_linear_projection(forest, subset = !is.na(X[, 1]))

  expect_equal(which.max(varimp), 1)
})
