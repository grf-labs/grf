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
  
  expect_true(min(var.big.group) > 0.025)
  expect_true(min(var.sm.group) > 0.025)
  expect_true(mean(var.big.group)/mean(var.sm.group) > 0.8)
})