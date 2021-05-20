library(grf)

set.seed(1234)

extract_samples <- function(tree) {

  # Keep only leaf nodes
  leaf_nodes <- Filter(f = function(x) x$is_leaf, tree$nodes)

  # Leaf nodes' 'samples' are estimation samples
  estimation_sample <- unlist(Map(f = function(x) x$samples, leaf_nodes))

  # Split = Drawn - Samples
  split_sample <- base::setdiff(tree$drawn_samples, estimation_sample)

  return(list(
    estimation_sample = estimation_sample,
    split_sample = split_sample
  ))
}

test_that("changing honest.fraction behaves as expected", {
  sample_fraction_1 <- 0.5
  honesty_fraction_1 <- 0.25

  sample_fraction_2 <- 0.25
  honesty_fraction_2 <- 0.1

  sample_fraction_3 <- 0.25
  honesty_fraction_3 <- 0.9

  n <- 16
  k <- 10
  X <- matrix(runif(n * k), nrow = n, ncol = k)
  Y <- runif(n)
  forest_1 <- grf::regression_forest(X, Y,
    sample.fraction = sample_fraction_1,
    honesty = TRUE, honesty.fraction = honesty_fraction_1
  )
  samples <- extract_samples(get_tree(forest_1, 1))

  expect_equal(length(samples$split_sample), n * sample_fraction_1 * honesty_fraction_1)
  expect_equal(length(samples$estimation_sample), n * sample_fraction_1 * (1 - honesty_fraction_1))
  # Checking for the runtime_error:
  # "The honesty fraction is too close to 1 or 0, as no observations will be sampled."
  expect_error(
    grf::regression_forest(X, Y,
      sample.fraction = sample_fraction_2,
      honesty = TRUE, honesty.fraction = honesty_fraction_2
    ),
    class = "std::runtime_error"
  )
  expect_error(
    grf::regression_forest(X, Y,
      sample.fraction = sample_fraction_3,
      honesty = TRUE, honesty.fraction = honesty_fraction_3
    ),
    class = "std::runtime_error"
  )
})

test_that("regression variance estimates are positive", {
  p <- 6
  n <- 1000

  ticks <- 101
  X.test <- matrix(0, ticks, p)
  xvals <- seq(-1, 1, length.out = ticks)
  X.test[, 1] <- xvals
  truth <- xvals > 0

  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- (X[, 1] > 0) + 2 * rnorm(n)

  forest <- regression_forest(X, Y, num.trees = 1000, ci.group.size = 4)
  preds.oob <- predict(forest, estimate.variance = TRUE)
  preds <- predict(forest, X.test, estimate.variance = TRUE)

  expect_true(all(preds$variance.estimate > 0))
  expect_true(all(preds.oob$variance.estimate > 0))

  error <- preds$predictions - truth
  expect_lt(mean(error^2), 0.2)

  truth.oob <- (X[, 1] > 0)
  error.oob <- preds.oob$predictions - truth.oob
  expect_lt(mean(error.oob^2), 0.2)

  Z.oob <- error.oob / sqrt(preds.oob$variance.estimate)
  expect_lt(mean(abs(Z.oob) > 1), 0.5)
})

test_that("Converting from a sparse data representation produces the same predictions", {
  dim <- 20
  X <- diag(rnorm(dim), dim)
  sparse.X <- Matrix::Matrix(X, sparse = TRUE)
  Y <- 1000 * (X[, 1]) + rnorm(dim)

  forest <- regression_forest(X, Y, mtry = dim, seed = 10)
  preds <- predict(forest, estimate.variance = TRUE)

  sparse.forest <- regression_forest(as.matrix(sparse.X), Y, mtry = dim, seed = 10)
  sparse.preds <- predict(sparse.forest, estimate.variance = TRUE)

  expect_equal(preds$predictions, sparse.preds$predictions)
})

test_that("OOB predictions contain debiased error estimates", {
  p <- 6
  n <- 10

  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- (X[, 1] > 0) + 2 * rnorm(n)

  forest <- regression_forest(X, Y, num.trees = 1000, ci.group.size = 4)
  preds.oob <- predict(forest)

  expect_equal(n, length(preds.oob$debiased.error))
})

test_that("regression forests with a positive imbalance.penalty have reasonable tree depths", {
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  Y <- 1000 * (X[, 1]) + rnorm(n)

  forest <- regression_forest(X, Y, imbalance.penalty = 0.001)
  split.freq <- split_frequencies(forest)
  expect_gt(sum(split.freq[4, ]), 0)
})

test_that("regression forests with a very small imbalance.penalty behave similarly to unpenalized forests.", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + 0.1 * rnorm(n)

  forest <- regression_forest(X, Y, imbalance.penalty = 0.0)
  forest.large.penalty <- regression_forest(X, Y, imbalance.penalty = 100.0)
  forest.small.penalty <- regression_forest(X, Y, imbalance.penalty = 1e-7)

  diff.large.penalty <- abs(forest.large.penalty$debiased.error - forest$debiased.error)
  diff.small.penalty <- abs(forest.small.penalty$debiased.error - forest$debiased.error)
  expect_lt(mean(diff.small.penalty), 0.10 * mean(diff.large.penalty))
})

test_that("variance estimates are positive [with sample weights]", {
  n <- 1000
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- abs(X[, 1]) + 0.1 * rnorm(n)
  e <- 1 / (1 + exp(-3 * X[, 1]))
  sample.weights <- 1 / e

  forest.weighted <- regression_forest(X, Y, sample.weights = sample.weights)
  mu.forest <- predict(forest.weighted, X, estimate.variance = TRUE)
  expect_true(all(mu.forest$variance.estimates > 0))
})

test_that("predictions are invariant to scaling of the sample weights.", {
  n <- 1000
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- abs(X[, 1]) + 0.1 * rnorm(n)
  e <- 1 / (1 + exp(-3 * X[, 1]))
  sample.weights <- 1 / e

  # The multiple is a power of 2 to avoid rounding errors allowing for exact comparison
  # between two forest with the same seed.
  forest.1 <- regression_forest(X, Y, sample.weights = sample.weights, seed = 1)
  forest.2 <- regression_forest(X, Y, sample.weights = 64 * sample.weights, seed = 1)
  pred.1 <- predict(forest.1, estimate.variance = TRUE)
  pred.2 <- predict(forest.2, estimate.variance = TRUE)

  expect_equal(pred.1$predictions, pred.2$predictions, tolerance = 1e-10)
  expect_equal(pred.1$variance.estimates, pred.2$variance.estimates, tolerance = 1e-10)
  expect_equal(pred.1$debiased.error, pred.2$debiased.error, tolerance = 1e-10)
})

test_that("sample weighted regression forest works as expected", {
  n <- 2000
  p <- 5
  obs.prob <- 1 / 20
  Y0 <- rbinom(n, 1, obs.prob / (1 + obs.prob))
  Y <- Y0 + rnorm(n) * 0.01
  X <- matrix(rnorm(n * p), n, p)
  sample.weights <- 1 + Y0 * (1 / obs.prob - 1)
  rf <- regression_forest(X, Y, sample.weights = sample.weights, num.trees = 250)

  expect_equal(mean(predict(rf)$predictions), weighted.mean(Y, sample.weights), tolerance = 0.05)
})

test_that("sample weighted regression forest is estimated with kernel weights `forest.weights * sample.weights`", {
  n <- 2000
  p <- 5
  obs.prob <- 1 / 20
  Y0 <- rbinom(n, 1, obs.prob / (1 + obs.prob))
  Y <- Y0 + rnorm(n) * 0.01
  X <- matrix(rnorm(n * p), n, p)
  sample.weights <- 1 + Y0 * (1 / obs.prob - 1)
  rf <- regression_forest(X, Y, sample.weights = sample.weights, num.trees = 250)

  x1 <- X[1, , drop = F]
  theta1 <- predict(rf, x1)$predictions
  alpha1 <- get_forest_weights(rf, x1)[1, ]
  theta1.lm <- lm(Y ~ 1, weights = alpha1 * sample.weights)

  expect_equal(theta1, theta1.lm$coefficients[[1]], tolerance = 1e-10)
})

test_that("sample weighting in the training of a regression forest improves its sample-weighted MSE.", {
  mse.ratio <- summary(replicate(4, {
    n <- 1000
    p <- 2
    X <- matrix(rnorm(n * p), n, p)
    Y <- abs(X[, 1]) + 0.1 * rnorm(n)
    e <- 1 / (1 + exp(-3 * X[, 1]))
    sample.weights <- 1 / e

    forest <- regression_forest(X, Y, num.trees = 500)
    forest.weighted <- regression_forest(X, Y, sample.weights = sample.weights, num.trees = 500)
    weighted.mse.forest <- sum(sample.weights * (forest$predictions - Y)^2)
    weighted.mse.forest.weighted <- sum(sample.weights * (forest.weighted$predictions - Y)^2)
    weighted.mse.forest.weighted / weighted.mse.forest
  }))
  expect_lt(mse.ratio[["1st Qu."]], 0.8)
})

test_that("inverse propensity weighting in the training of a regression forest with missing data improves
           its complete-data MSE.", {
  mse.ratio <- summary(replicate(4, {
    n <- 1000
    p <- 2
    X <- matrix(rnorm(n * p), n, p)
    Y <- abs(X[, 1]) + 0.1 * rnorm(n)
    e <- 1 / (1 + exp(-3 * X[, 1]))
    w <- runif(n) <= e
    sample.weights <- 1 / e[w]

    forest <- regression_forest(X[w, ], Y[w], num.trees = 500)
    forest.weighted <- regression_forest(X[w, ], Y[w], sample.weights = sample.weights, num.trees = 500)
    ipw.mse.forest <- sum((predict(forest, X) - Y)^2)
    ipw.mse.forest.weighted <- sum((predict(forest.weighted, X) - Y)^2)
    ipw.mse.forest.weighted / ipw.mse.forest
  }))
  # Due to the influence of rare extreme weights we check the following more conservative bound.
  expect_lt(mse.ratio[["1st Qu."]], 0.85)
})

test_that("sample weighted regression forest gives correct coverage", {
  n <- 2500
  p <- 5
  pA <- 0.2
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, pA)
  gamma <- A / pA + (1 - A) / (1 - pA)
  Y.true <- 1 / (1 + exp(-3 * X[,1]))
  Y <- 2 * Y.true * A + rnorm(n)

  rf <- regression_forest(X, Y, sample.weights = gamma)
  preds <- predict(rf, estimate.variance = TRUE)
  zstat <- (preds$predictions - Y.true) / sqrt(preds$variance.estimates)
  coverage <- mean(abs(zstat) < 1.96)
  mse <- mean((preds$predictions - Y.true)^2)

  rf.noweight <- regression_forest(X, Y)
  preds.noweight <- predict(rf.noweight, estimate.variance = TRUE)
  zstat.noweight <- (preds.noweight$predictions - Y.true) / sqrt(preds.noweight$variance.estimates)
  coverage.noweight <- mean(abs(zstat.noweight) < 1.96)
  mse.noweight <- mean((preds.noweight$predictions - Y.true)^2)

  expect_lt(mse / mse.noweight, 0.3)
  expect_gt(coverage, 0.8)
  expect_lt(coverage.noweight, 0.55)
})

test_that("sample weighted regression forest is identical to replicating samples", {
  # To make these forests comparable sample.fraction has to be 1 to draw the same samples
  # and min.node.size 1 for the split stopping condition to be the same.
  n <- 500
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)
  to.duplicate <- sample(1:n, 100)
  XX <- rbind(X, X[to.duplicate, ])
  YY <- c(Y, Y[to.duplicate])
  sample.weights <- rep(1, n)
  sample.weights[to.duplicate] <- 2

  rf.weighted <- regression_forest(X, Y, sample.weights = sample.weights,
                                   num.trees = 500,
                                   sample.fraction = 1,
                                   min.node.size = 1,
                                   honesty = FALSE,
                                   ci.group.size = 1,
                                   seed = 123)

  rf.duplicated.data <- regression_forest(XX, YY,
                                          num.trees = 500,
                                          sample.fraction = 1,
                                          min.node.size = 1,
                                          honesty = FALSE,
                                          ci.group.size = 1,
                                          seed = 123)
  expect_equal(predict(rf.weighted, XX)$predictions, predict(rf.duplicated.data, XX)$predictions)
})

test_that("a non-pruned honest regression forest has lower MSE than a pruned honest regression forests
          (on small data)", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- abs(X[, 1]) + 0.1 * rnorm(n)

  f1 <- regression_forest(X, Y, honesty = TRUE, honesty.fraction = 0.9, honesty.prune.leaves = TRUE)
  f2 <- regression_forest(X, Y, honesty = TRUE, honesty.fraction = 0.9, honesty.prune.leaves = FALSE)

  mse.pruned <- mean((predict(f1)$predictions - Y)^2)
  mse.notpruned <- mean((predict(f2)$predictions - Y)^2)

  # Upper bound of 65 % is based on 10 000 repetitions of the above DGP
  expect_lt(mse.notpruned / mse.pruned, 0.65)
})

test_that("regression_forest works as expected with missing values", {
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)

  nmissing <- sample(c(1, 150, 500, 1000), size = 1)
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

  rf.mia <- regression_forest(X.mia, Y, seed = 123)
  rf <- regression_forest(X, Y, seed = 123)

  mse.oob.diff <- mean((predict(rf.mia)$predictions - predict(rf)$predictions)^2)
  mse.diff <- mean((predict(rf.mia, X.mia)$predictions - predict(rf, X)$predictions)^2)
  mse.test.diff <- mean((predict(rf.mia, X.mia.test)$predictions - predict(rf, X.test)$predictions)^2)

  diff.mse.oob <- mean((predict(rf.mia)$predictions - Y)^2) - mean((predict(rf)$predictions - Y)^2)
  diff.mse <- mean((predict(rf.mia, X.mia)$predictions - Y)^2) - mean((predict(rf, X)$predictions - Y)^2)
  diff.mse.test <- mean((predict(rf.mia, X.mia.test)$predictions - Y)^2) - mean((predict(rf, X.test)$predictions - Y)^2)

  expect_equal(mse.oob.diff, 0, tolerance = 0.001)
  expect_equal(mse.diff, 0, tolerance = 0.001)
  expect_equal(mse.test.diff, 0, tolerance = 0.001)

  expect_equal(diff.mse.oob, 0, tolerance = 0.01)
  expect_equal(diff.mse, 0, tolerance = 0.01)
  expect_equal(diff.mse.test, 0, tolerance = 0.01)

  # All NaNs
  X[, ] <- NaN

  Xl <- X
  Xr <- X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)

  rf.mia <- regression_forest(X.mia, Y, seed = 123)
  rf <- regression_forest(X, Y, seed = 123)
  mse.oob.diff.allnan <- mean((predict(rf.mia)$predictions - predict(rf)$predictions)^2)
  expect_equal(mse.oob.diff.allnan, 0, tolerance = 0.0001)
})
