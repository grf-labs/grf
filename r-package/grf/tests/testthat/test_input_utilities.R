library(grf)

test_that("forests handle matrix and data.frame inputs equally", {
  
  n = 100
  p = 2
  Xm = matrix(rnorm(n*p), n, p)
  Xd = as.data.frame(Xm)
  Y = rnorm(n)
  W = rbinom(0.5, size=1, n=n)
  Z = runif(n)
  
  # Matrix input
  rfm = regression_forest(Xm, Y, num.trees=5, seed=1234)
  cfm = causal_forest(Xm, Y, W, num.trees=5, seed=1234)
  ifm = instrumental_forest(Xm, Y, W, Z, num.trees=5, seed=1234)
  
  # Data.frame input
  rfd = regression_forest(Xd, Y, num.trees=5, seed=1234)
  cfd = causal_forest(Xd, Y, W, num.trees=5, seed=1234)
  ifd = instrumental_forest(Xd, Y, W, Z, num.trees=5, seed=1234)
  
  # Check that output is the same
  expect_equal(predict(rfd)$predictions, predict(rfm)$predictions)
  expect_equal(predict(cfd)$predictions, predict(cfm)$predictions)
  expect_equal(predict(ifd)$predictions, predict(ifm)$predictions)

})
