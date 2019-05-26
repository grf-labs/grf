library(grf)
library(testthat)

test_that("forests take both matrix and data.frame inputs", {
  
  n = 100
  p = 2
  Xm = matrix(rnorm(n*p), n, p)
  Xd = as.data.frame(Xm)
  Y = rnorm(n)
  W = rbinom(0.5, size=1, n=n)
  Z = runif(n)
  
  # Matrix input
  rf = regression_forest(Xm, Y, num.trees=5)
  cf = causal_forest(Xm, Y, W, num.trees=5)
  zf = instrumental_forest(Xm, Y, W, Z, num.trees=5)
  
  # Data.frame input
  rf = regression_forest(Xd, Y, num.trees=5)
  cf = causal_forest(Xd, Y, W, num.trees=5)
  zf = instrumental_forest(Xd, Y, W, Z, num.trees=5)
  
})
