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
  expect_true(mean(error^2) < 0.2)

  truth.oob <- (X[, 1] > 0)
  error.oob <- preds.oob$predictions - truth.oob
  expect_true(mean(error.oob^2) < 0.2)

  Z.oob <- error.oob / sqrt(preds.oob$variance.estimate)
  expect_true(mean(abs(Z.oob) > 1) < 0.5)
})

test_that("using a sparse data representation produces the same predictions", {
  dim <- 20
  X <- diag(rnorm(dim), dim)
  sparse.X <- as(X, "dgCMatrix")
  Y <- 1000 * (X[, 1]) + rnorm(dim)

  forest <- regression_forest(X, Y, mtry = dim, seed = 10)
  preds <- predict(forest, estimate.variance = TRUE)

  sparse.forest <- regression_forest(sparse.X, Y, mtry = dim, seed = 10)
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
  expect_true(sum(split.freq[4, ]) > 0)
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
  expect_true(mean(diff.small.penalty) < 0.10 * mean(diff.large.penalty))
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

test_that("debiased errors are smaller than raw errors [with sample weights]", {
  n <- 1000
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- abs(X[, 1]) + 0.1 * rnorm(n)
  e <- 1 / (1 + exp(-3 * X[, 1]))
  sample.weights <- 1 / e

  forest <- regression_forest(X, Y, sample.weights = sample.weights)
  preds <- predict(forest)
  expect_true(all(preds$debiased.error^2 < preds$error^2))
})

test_that("predictions are invariant to scaling of the sample weights.", {
  n <- 1000
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- abs(X[, 1]) + 0.1 * rnorm(n)
  e <- 1 / (1 + exp(-3 * X[, 1]))
  sample.weights <- 1 / e

  forest.1 <- regression_forest(X, Y, sample.weights = sample.weights)
  forest.2 <- regression_forest(X, Y, sample.weights = 1e-6 * sample.weights)
  expect_true(max(abs(forest.1$predictions - forest.2$predictions)) < .1)
  # forests are built with different random seeds, hence possibly poor agreement
})

test_that("sample weighting in the training of a regression forest improves its sample-weighted MSE.", {
  n <- 1000
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- abs(X[, 1]) + 0.1 * rnorm(n)
  e <- 1 / (1 + exp(-3 * X[, 1]))
  sample.weights <- 1 / e

  forest <- regression_forest(X, Y)
  forest.weighted <- regression_forest(X, Y, sample.weights = sample.weights)
  weighted.mse.forest <- sum(sample.weights * (forest$predictions - Y)^2)
  weighted.mse.forest.weighted <- sum(sample.weights * (forest.weighted$predictions - Y)^2)
  expect_true(weighted.mse.forest.weighted < weighted.mse.forest)
})

test_that("inverse propensity weighting in the training of a regression forest with missing data improves
           its complete-data MSE.", {
  n <- 1000
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- abs(X[, 1]) + 0.1 * rnorm(n)
  e <- 1 / (1 + exp(-3 * X[, 1]))
  w <- runif(n) <= e
  sample.weights <- 1 / e[w]

  forest <- regression_forest(X[w, ], Y[w])
  forest.weighted <- regression_forest(X[w, ], Y[w], sample.weights = sample.weights)
  ipw.mse.forest <- sum((predict(forest, X) - Y)^2)
  ipw.mse.forest.weighted <- sum((predict(forest.weighted, X) - Y)^2)
  expect_true(ipw.mse.forest.weighted < ipw.mse.forest)
})

test_that("sample weighting is identical to replicating samples", {
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
  diff.abs <- abs(predict(rf.weighted, XX)$predictions - predict(rf.duplicated.data, XX)$predictions)

  expect_equal(all(diff.abs == 0), TRUE)
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
  expect_true(mse.notpruned < 0.65 * mse.pruned)
})
