library(grf)

set.seed(123)

test_that("using big ci.group.size doesn't result in tiny CIs", {
  p = 4
  n = 200
  
  X = matrix(2 * runif(n * p) - 1, n, p)
  Y = (X[,1] > 0) + rnorm(n)
  X.test = matrix(2 * runif(2000 * p) - 1, 2000, p)
  
  forest.big.group = regression_forest(X, Y, num.trees = 2000, ci.group.size = 200)
  pred.big.group = predict(forest.big.group, X.test, estimate.variance=TRUE)
  var.big.group = sqrt(pred.big.group$variance.estimates)
  
  forest.sm.group = regression_forest(X, Y, num.trees = 2000, ci.group.size = 10)
  pred.sm.group = predict(forest.sm.group, X.test, estimate.variance=TRUE)
  var.sm.group = sqrt(pred.sm.group$variance.estimates)

  expect_true(min(var.big.group) > 0.01)
  expect_true(min(var.sm.group) > 0.01)
  expect_true(mean(var.big.group)/mean(var.sm.group) > 0.8)
})


test_that("regression CIs are reasonable", {
  p = 6
  n = 1000
  X = matrix(2 * runif(n * p) - 1, n, p)
  Y = (X[,1] > 0) + 2 * rnorm(n)
  forest = regression_forest(X, Y)
  preds.oob = predict(forest, estimate.variance = TRUE)
  error.standardized = (preds.oob$predictions - (X[,1] > 0)) / sqrt(preds.oob$variance.estimates)
  expect_true(mean(abs(error.standardized) > qnorm(0.975)) <= 0.15)
})

test_that("instrumental CIs are reasonable", {
  n = 1000
  n.test = 1000
  p = 4
  alpha.mu = 3
  alpha.tau = 1
  k.tau = 2
  k.mu = 2
  
  X = matrix(rnorm(n * p), n, p)
  eps = 2 * rnorm(n)
  Z = rbinom(n, 1, 2/3)
  filter = rbinom(n, 1, 1/(1 + exp(-eps)))
  W = Z * filter
  tau = alpha.tau * apply(X[,1:k.tau], 1, function(xx) sum(pmax(0, xx)))
  mu = alpha.mu *  apply(X[,k.tau + 1:k.mu], 1, function(xx) sum(pmax(0, xx)))
  Y = (2 * W - 1) / 2 * tau + mu + eps
  
  X.test = matrix(rnorm(n.test * p), n.test, p)
  tau.true = alpha.tau *  apply(X.test[,1:k.tau], 1, function(xx) sum(pmax(0, xx)))
  
  forest = instrumental_forest(X, Y, W, Z, precompute.nuisance = TRUE)
  tau.hat = predict(forest, newdata = X.test, estimate.variance = TRUE)
  error.standardized = (tau.hat$predictions - tau.true) / sqrt(tau.hat$variance.estimates)
  expect_true(mean(abs(error.standardized) > qnorm(0.975)) <= 0.18)
})

test_that("instrumental CIs are invariant to scaling Z", {
  n = 2000
  p = 5
  n.test = 2000
  X = matrix(rnorm(n*p), n, p)
  W = rnorm(n)
  Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
  
  X.test = matrix(rnorm(n*p), n, p)
  tau.true = pmax(X.test[,1], 0)
  
  forest = causal_forest(X, Y, W, precompute.nuisance = FALSE)
  tau.hat = predict(forest, newdata = X.test, estimate.variance = TRUE)
  error.standardized = (tau.hat$predictions - tau.true) / sqrt(tau.hat$variance.estimates)
  expect_true(mean(abs(error.standardized) > qnorm(0.975)) <= 0.15)
  expect_true(mean(abs(error.standardized) > qnorm(0.975)) >= 0.005)
  
  Z = 0.00000001 * W
  forest.iv = instrumental_forest(X, Y, W, Z, precompute.nuisance = FALSE)
  tau.hat.iv = predict(forest.iv, newdata = X.test, estimate.variance = TRUE)
  error.standardized.iv = (tau.hat.iv$predictions - tau.true) / sqrt(tau.hat.iv$variance.estimates)
  expect_true(mean(abs(error.standardized.iv) > qnorm(0.975)) <= 0.15)
  expect_true(mean(abs(error.standardized.iv) > qnorm(0.975)) >= 0.005)
  
  kst = ks.test(error.standardized, error.standardized.iv)
  expect_true(kst$statistic <= 0.05)
})
