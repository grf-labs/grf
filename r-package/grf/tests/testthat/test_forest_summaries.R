library(grf)

set.seed(1234)

test_that("causal forest calibration is reasonable", {
    n = 800; p = 4
    X = matrix(rnorm(n*p), n, p)
    W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
    Y = pmax(X[,1], 0) * (W - 0.75) + rnorm(n)

    cf = causal_forest(seed=1000, X, Y, W,
                       W.hat = 0.25 + 0.5 * (X[,1] > 0),
                       Y.hat = pmax(X[,1], 0) * (0.5 * (X[,1] > 0) - 0.5),
                       num.trees = 500)
    tc = test_calibration(cf)
    
    expect_true(abs(tc[1,1] - 1) <= 0.4)
    expect_true(abs(tc[2,1] - 1) <= 0.4)
})

test_that("causal forest calibration is reasonable with no average effect", {
  n = 800; p = 4
  X = matrix(rnorm(n*p), n, p)
  W = rnorm(n, 1, 0.5)
  Y = sign(X[,1]) * (W - 0.5) + rnorm(n)
  
  cf = causal_forest(seed=1000, X, Y, W,
                     W.hat = 0.5,
                     Y.hat = 0,
                     num.trees = 500)
  tc = test_calibration(cf)
  
  expect_true(abs(tc[1,3]) <= 3)
  expect_true(abs(tc[2,1] - 1) <= 0.3)
})

test_that("causal forest calibration is reasonable with no heterogeneous effect", {
  n = 800; p = 4
  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
  Y = pmax(X[,2], 0) + W + rnorm(n)

  cf = causal_forest(seed=1000, X, Y, W,
                     W.hat = 0.25 + 0.5 * (X[,1] > 0),
                     Y.hat = 0.25 + 0.5 * (X[,1] > 0) + pmax(X[,2], 0),
                     num.trees = 500)
  tc = test_calibration(cf)

  expect_true(abs(tc[1,1] - 1) <= 0.3)
  expect_true(abs(tc[2,3]) <= 4)
})

test_that("causal forest calibration is reasonable with no heterogeneous effect with sample weights and clusters", {
  p = 4
  K = 100
  cluster.sizes = pmax(20, round(40+3*rt(K, df=3)))
  n = sum(cluster.sizes)
  clust = rep(1:K, cluster.sizes)

  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
  big = as.numeric(cluster.sizes[clust] >= median(cluster.sizes))
  Y = pmax(X[,2], 0)*big + W + rnorm(n)

  e.cc = 1/( 1+exp(-2*X[,1])  )
  cc = as.logical(rbinom(n, 1, e.cc))
  sample.weights = 1/e.cc

  cf = causal_forest(seed=1000, X[cc,], Y[cc], W[cc],
                     W.hat = 0.25 + 0.5 * (X[cc,1] > 0),
                     Y.hat = 0.25 + 0.5 * (X[cc,1] > 0) + pmax(X[cc,2], 0) * big[cc],
                     sample.weights = sample.weights[cc],
                     clusters=clust[cc],
                     num.trees = 500)
  tc = test_calibration(cf)

  expect_true(abs(tc[1,1] - 1) <= 0.3)
  expect_true(abs(tc[2,3]) <= 4)
})

test_that("regression forest calibration is reasonable", {
  n = 100; p = 4
  X = matrix(rnorm(n*p), n, p)
  Y = 5 + 5 * sign(X[,1]) + rnorm(n)
  
  rf = regression_forest(seed=1000, X, Y)
  tc = test_calibration(rf)
  
  expect_true(abs(tc[1,1] - 1) <= 0.1)
  expect_true(abs(tc[2,1] - 1) <= 0.2)
})

test_that("regression forest calibration is reasonable with no heterogeneous effect", {
  n = 100; p = 4
  X = matrix(rnorm(n*p), n, p)
  Y = 5 + rnorm(n)
  
  rf = regression_forest(seed=1000, X, Y)
  tc = test_calibration(rf)
  
  expect_true(abs(tc[1,1] - 1) <= 0.1)
  expect_true(abs(tc[2,3]) <= 4)
})

test_that("causal forest calibration works with clusters", {
  n = 100; p = 4
  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
  Y = pmax(X[,2], 0) + W + rnorm(n)
  
  cf = causal_forest(seed=1000, X, Y, W,
                     W.hat = 0.25 + 0.5 * (X[,1] > 0),
                     Y.hat = 0.25 + 0.5 * (X[,1] > 0) + pmax(X[,2], 0),
                     num.trees = 100)
  tc = test_calibration(cf)
  
  cf.clust = cf
  cf.clust$W.orig = c(cf$W.orig[1:(n/2)], rep(cf$W.orig[n/2 + 1:(n/2)], 10))
  cf.clust$Y.orig = c(cf$Y.orig[1:(n/2)], rep(cf$Y.orig[n/2 + 1:(n/2)], 10))
  cf.clust$W.hat = c(cf$W.hat[1:(n/2)], rep(cf$W.hat[n/2 + 1:(n/2)], 10))
  cf.clust$Y.hat = c(cf$Y.hat[1:(n/2)], rep(cf$Y.hat[n/2 + 1:(n/2)], 10))
  cf.clust$predictions = c(cf$predictions[1:(n/2)], rep(cf$predictions[n/2 + 1:(n/2)], 10))
  cf.clust$debiased.error = c(cf$debiased.error[1:(n/2)], rep(cf$debiased.error[n/2 + 1:(n/2)], 10))
  cf.clust$excess.error = c(cf$excess.error[1:(n/2)], rep(cf$excess.error[n/2 + 1:(n/2)], 10))
  cf.clust$clusters = c(1:(n/2), rep(n/2 + 1:(n/2), 10))
  tc.clust = test_calibration(cf.clust)
  
  expect_equal(tc[,1:3], tc.clust[,1:3])
})
