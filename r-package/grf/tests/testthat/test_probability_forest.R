library(grf)

test_that("probability forest works as expected", {
  p <- 5
  n <- 200
  X <- matrix(rnorm(n * p), n, p)
  Y <- as.factor(sample(c(0, 1, 5, 8), n, T))
  class.names <- levels(Y)

  prf <- probability_forest(X, Y, num.trees = 50)
  pred.oob <- predict(prf, estimate.variance = TRUE)
  pred <- predict(prf, X, estimate.variance = TRUE)

  expect_equal(colnames(pred$predictions), class.names)
  expect_equal(rowSums(pred$predictions), rep(1, n), tolerance = 1e-10)
  expect_true(all(pred$predictions >= 0))
  expect_true(all(pred$variance.estimates > 0))

  expect_equal(colnames(pred.oob$predictions), class.names)
  expect_equal(rowSums(pred.oob$predictions), rep(1, n), tolerance = 1e-10)

  Y <- sample(c("A", "B", "C"), n, T)
  Y <- as.factor(Y)
  class.names <- levels(Y)
  prf <- probability_forest(X, Y, num.trees = 50)
  pred.oob <- predict(prf)
  expect_equal(colnames(pred.oob$predictions), class.names)

  p <- 5
  n <- 200
  X <- matrix(rnorm(n*p), n, p)
  prob <- 1 / (1 + exp(-X[, 1] - X[, 2]))
  Y <- rbinom(n, 1, prob)
  pf <- probability_forest(X, as.factor(Y), seed = 42, num.trees = 100)
  pf.flipped <- probability_forest(X, as.factor(1 - Y), seed = 42, num.trees = 100)

  expect_equal(unname(predict(pf)$predictions), unname(predict(pf.flipped)$predictions[, c(2, 1)]))
  expect_equal(unname(predict(pf, X)$predictions), unname(predict(pf.flipped, X)$predictions[, c(2, 1)]))
})

test_that("probability forest is well-calibrated", {
  p <- 5
  n <- 500
  X <- matrix(runif(n * p), n, p)
  X <- sweep(X, 1, rowSums(X), "/")
  prob.true <- X[, 1:4]
  Y <- sapply(1:n, function(i) sample(c("a", "b", "c", "d"), 1, prob = prob.true[i, ]))
  Y <- as.factor(Y)

  prf <- probability_forest(X, Y, num.trees = 500)
  p.hat <- predict(prf, estimate.variance = TRUE)
  expect_lt(mean(rowMeans((p.hat$predictions - prob.true)^2)), 0.01)

  z.score <- abs(p.hat$predictions - prob.true) / sqrt(p.hat$variance.estimates)
  coverage <- colMeans(z.score <= 1.96)
  expect_gt(coverage[1], 0.7)
  expect_gt(coverage[2], 0.7)
  expect_gt(coverage[3], 0.7)
  expect_gt(coverage[4], 0.7)
})

test_that("sample weighted probability forest is invariant to scaling", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
  Y <- as.factor(Y)
  weights <- runif(n)

  # Predictions are invariant to sample weight scaling
  pf <- probability_forest(X, Y, num.trees = 100, sample.weights = weights, seed = 1)
  pf2 <- probability_forest(X, Y, num.trees = 100, sample.weights = weights * 2, seed = 1)
  p1 <- predict(pf)
  p2 <- predict(pf2)

  expect_equal(p1$predictions, p2$predictions)

  # Variance estimates are invariant to sample weight scaling
  v1 <- predict(pf, estimate.variance = TRUE)
  v2 <- predict(pf2, estimate.variance = TRUE)

  expect_equal(v1$variance.estimates, v2$variance.estimates)
})

test_that("sample weighted probability forest is identical to replicating samples", {
  # To make these forests comparable sample.fraction has to be 1 to draw the same samples
  # and min.node.size 1 for the split stopping condition to be the same.
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-5 * X[, 1]))
  Y <- rbinom(n, 1, e)
  to.duplicate <- sample(1:n, 100)
  XX <- rbind(X, X[to.duplicate, ])
  YY <- c(Y, Y[to.duplicate])
  sample.weights <- rep(1, n)
  sample.weights[to.duplicate] <- 2

  pf.weighted <- probability_forest(X, as.factor(Y),
                                    sample.weights = sample.weights,
                                    num.trees = 250,
                                    sample.fraction = 1,
                                    min.node.size = 1,
                                    honesty = FALSE,
                                    ci.group.size = 1,
                                    seed = 42)

  pf.duplicated <- probability_forest(XX, as.factor(YY),
                                      num.trees = 250,
                                      sample.fraction = 1,
                                      min.node.size = 1,
                                      honesty = FALSE,
                                      ci.group.size = 1,
                                      seed = 42)

  expect_equal(predict(pf.weighted, XX)$predictions, predict(pf.duplicated, XX)$predictions)
})

test_that("sample weighted probability forest improves complete-data MSE", {
  mean <- mean(replicate(4, {
      n <- 500
      p <- 5
      X <- matrix(rnorm(n * p), n, p)
      e <- 1 / (1 + exp(-5 * X[, 1]))
      Y <- as.factor(rbinom(n, 1, e))
      cc <- runif(n) < e
      sample.weights <- 1 / e[cc]
      pp.true <- cbind(1 - e, e)

      pf <- probability_forest(X[cc, ], Y[cc], num.trees = 500)
      pf.weighted <- probability_forest(X[cc, ], Y[cc], sample.weights = sample.weights, num.trees = 500)

      mse <- mean((predict(pf, X)$predictions - pp.true)^2)
      mse.weighted <- mean((predict(pf.weighted, X)$predictions - pp.true)^2)

      mse.weighted / mse
    }))

  expect_lt(mean, 0.8)
})

test_that("sample weighted probability forest works as expected", {
  n <- 2000
  p <- 5
  obs.prob <- 1 / 20
  Y <- rbinom(n, 1, obs.prob / (1 + obs.prob))
  X <- matrix(rnorm(n * p), n, p)
  sample.weights <- 1 + Y * (1 / obs.prob - 1)

  pf <- probability_forest(X, as.factor(Y), sample.weights = sample.weights, num.trees = 500)
  pp <- predict(pf)$predictions

  expect_equal(mean(pp[, "1"]), weighted.mean(Y, sample.weights), tolerance = 0.05)
})

test_that("probability forest with sample weights is kernel weighted correctly", {
  p <- 5
  n <- 500
  X <- matrix(rnorm(n*p), n, p)
  prob <- 1 / (1 + exp(-X[, 1] - X[, 2]))
  Y <- as.factor(rbinom(n, 1, prob))
  sample.weights <- sample(c(1, 10), n, TRUE)
  pf <- probability_forest(X, Y, sample.weights = sample.weights, num.trees = 250)

  x1 <- X[1, , drop = FALSE]
  theta1 <- predict(pf, x1)$predictions[,]

  alpha1 <- get_forest_weights(pf, x1)[1, ]
  y <- as.numeric(Y) - 1
  sums <- cbind((y == 0) * alpha1 * sample.weights, (y == 1) * alpha1 * sample.weights)
  theta1.weighted <- colSums(sums) / sum(sums)

  expect_equal(unname(theta1), theta1.weighted)
})
