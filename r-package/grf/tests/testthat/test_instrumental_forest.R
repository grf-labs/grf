library(grf)

set.seed(1234)

test_that("instrumental forests give reasonable estimates", {
  p = 6
  n = 2000
  n.test = 2000
  
  X = matrix(rnorm(n * p), n, p)
  
  eps = rnorm(n)
  Z = rbinom(n, 1, 2/3)
  filter = rbinom(n, 1, 1/(1 + exp(-1 * eps)))
  W = Z * filter
  
  tau = apply(X[,1:2], 1, function(xx) sum(pmax(0, xx)))
  mu = apply(X[,2 + 1:2], 1, function(xx) sum(pmax(0, xx)))
  
  Y = (2 * W - 1) / 2 * tau + mu + eps

  X.test = matrix(rnorm(n.test * p), n.test, p)  
  tau.true = apply(X.test[,1:2], 1, function(xx) sum(pmax(0, xx)))
  
  forest.iv = instrumental_forest(X, Y, W, Z, num.trees = 4000)
  
  preds.iv.oob = predict(forest.iv, estimate.variance=TRUE)
  preds.iv = predict(forest.iv, X.test, estimate.variance=TRUE)
  
  expect_true(all(preds.iv$variance.estimate > 0))
  expect_true(all(preds.iv.oob$variance.estimate > 0))
  
  error = preds.iv$predictions - tau.true
  expect_true(mean(error^2) < 0.4)
  
  error.oob = preds.iv.oob$predictions - tau
  expect_true(mean(error.oob^2) < 0.4)
  
  Z = error / sqrt(preds.iv$variance.estimate)
  Z.oob = error.oob / sqrt(preds.iv.oob$variance.estimate)
  
  expect_true(mean(abs(Z) > 1) < 0.5)
  expect_true(mean(abs(Z.oob) > 1) < 0.5)

})

test_that("instrumental forests with censoring and ipcc weights compare to forests given full data as you would expect:
 acceptable error, greater variance, comparable validity of CIs", {
  p = 6
  n = 2000
  n.test = 2000

  X = matrix(rnorm(n * p), n, p)

  eps = rnorm(n)
  Z = rbinom(n, 1, 2/3)
  filter = rbinom(n, 1, 1/(1 + exp(-1 * eps)))
  W = Z * filter

  tau = apply(X[,1:2], 1, function(xx) sum(pmax(0, xx)))
  mu = apply(X[,2 + 1:2], 1, function(xx) sum(pmax(0, xx)))

  Y = (2 * W - 1) / 2 * tau + mu + eps

  e.cc = 1/(1+exp(-3*X[,1]))
  cc = as.logical(rbinom(n, 1, e.cc))
  sample.weights = 1/e.cc

  X.test = matrix(rnorm(n.test * p), n.test, p)
  tau.true = apply(X.test[,1:2], 1, function(xx) sum(pmax(0, xx)))

  forest.iv.full =   forest.iv.ipcc = instrumental_forest(X, Y, W, Z, num.trees = 4000)
  forest.iv.ipcc = instrumental_forest(X[cc,], Y[cc], W[cc], Z[cc],
    sample.weights = sample.weights, num.trees = 4000)

  preds.iv.full = predict(forest.iv.full, X.test, estimate.variance=TRUE)
  preds.iv.ipcc = predict(forest.iv.ipcc, estimate.variance=TRUE)

  error.full = preds.iv.full$predictions - tau.true
  error.ipcc = preds.iv.ipcc$predictions - tau.true
  expect_true(mean(error.full^2) < 0.4)
  expect_true(mean(error.ipcc^2) < 0.4)

  expect_true(mean(preds.iv.full$variance.estimate / preds.iv.ipcc$variance.estimate) <= 1)

  Z.full = error.full / sqrt(preds.iv.full$variance.estimate)
  Z.ipcc = error.ipcc / sqrt(preds.iv.ipcc$variance.estimate)

  expect_true(mean(abs(Z.full) > 1) < 0.5)
  expect_true(mean(abs(Z.ipcc) > 1) < 0.5)
})