library(grf)

test_that("probability forest works as expected", {
  p <- 5
  n <- 200
  X <- matrix(rnorm(n * p), n, p)
  Y <- sample(c(0, 1, 5, 8), n, T)
  class.names <- sort(unique(Y))

  prf <- probability_forest(X, Y, num.trees = 50)
  pred.oob <- predict(prf, estimate.variance = TRUE)
  pred <- predict(prf, X, estimate.variance = TRUE)

  expect_true(all(colnames(pred$predictions) == class.names))
  expect_true(all(colnames(pred$variance.estimates) == class.names))
  expect_true(all(abs(rowSums(pred$predictions) - 1) < 1e-10))
  expect_true(all(pred$predictions >= 0))
  expect_true(all(pred$variance.estimates > 0))

  expect_true(all(colnames(pred.oob$predictions) == class.names))
  expect_true(all(colnames(pred.oob$variance.estimates) == class.names))
  expect_true(all(abs(rowSums(pred.oob$predictions) - 1) < 1e-10))

  Y <- sample(c("A", "B", "C"), n, T)
  class.names <- sort(unique(Y))
  prf <- probability_forest(X, Y, num.trees = 50)
  pred.oob <- predict(prf)
  expect_true(all(colnames(pred.oob$predictions) == class.names))
})

test_that("probability forest is well-calibrated", {
  p <- 5
  n <- 500
  X <- matrix(runif(n * p), n, p)
  X <- sweep(X, 1, rowSums(X), "/")
  prob.true <- X[, 1:4]
  Y <- sapply(1:n, function(i) sample(c("a", "b", "c", "d"), 1, prob = prob.true[i, ]))

  prf <- probability_forest(X, Y, num.trees = 500)
  p.hat <- predict(prf, estimate.variance = TRUE)

  expect_true(mean(rowMeans((p.hat$predictions - prob.true)^2)) < 0.01)

  ub <- p.hat$predictions + 2 * sqrt(p.hat$variance.estimates)
  lb <- p.hat$predictions - 2 * sqrt(p.hat$variance.estimates)
  covered <- lb < prob.true & prob.true < ub

  expect_true(all(colMeans(covered) > 0.7))
})
