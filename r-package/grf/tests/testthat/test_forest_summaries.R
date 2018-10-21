library(grf)

set.seed(1234)

test_that("causal forest calibration is reasonable", {
    n = 800; p = 4
    X = matrix(rnorm(n*p), n, p)
    W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
    Y = pmax(X[,1], 0) * (W - 0.75) + rnorm(n)

    cf = causal_forest(X, Y, W,
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
  
  cf = causal_forest(X, Y, W,
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
  
  cf = causal_forest(X, Y, W,
                     W.hat = 0.25 + 0.5 * (X[,1] > 0),
                     Y.hat = 0.25 + 0.5 * (X[,1] > 0) + pmax(X[,2], 0),
                     num.trees = 500)
  tc = test_calibration(cf)
  
  expect_true(abs(tc[1,1] - 1) <= 0.3)
  expect_true(abs(tc[2,3]) <= 4)
})

test_that("regression forest calibration is reasonable", {
  n = 100; p = 4
  X = matrix(rnorm(n*p), n, p)
  Y = 5 + 5 * sign(X[,1]) + rnorm(n)
  
  rf = regression_forest(X, Y)
  tc = test_calibration(rf)
  
  expect_true(abs(tc[1,1] - 1) <= 0.1)
  expect_true(abs(tc[2,1] - 1) <= 0.2)
})

test_that("regression forest calibration is reasonable with no heterogeneous effect", {
  n = 100; p = 4
  X = matrix(rnorm(n*p), n, p)
  Y = 5 + rnorm(n)
  
  rf = regression_forest(X, Y)
  tc = test_calibration(rf)
  
  expect_true(abs(tc[1,1] - 1) <= 0.1)
  expect_true(abs(tc[2,3]) <= 4)
})