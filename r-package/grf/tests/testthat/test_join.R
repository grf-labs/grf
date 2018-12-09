library(grf)

test_that("Concatenated regression forest attributes are sensible", {
  # Train regression forests
  n = 50; p = 2
  X = matrix(rnorm(n*p), n, p)
  Y = X[,1] * rnorm(n)
  W = X[,2] > 0
  r.forest1 = regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10)
  r.forest2 = regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 10)

  # Join the forests together. 
  big_rf = join_forests(list(r.forest1, r.forest2))

  # Result is also a regression_forest of the same class
  expect_true(big_rf$num.trees == (r.forest1$num.trees + r.forest2$num.trees))
  expect_true(is(big_rf, "grf"))
  expect_equal(class(r.forest1), class(big_rf))
})


test_that("Concatenated causal forest attributes are sensible", {
  # Train causal forests
  n = 50; p = 2
  X = matrix(rnorm(n*p), n, p)
  Y = X[,1] * rnorm(n)
  W = X[,2] > 0
  c.forest1 = causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 10)
  c.forest2 = causal_forest(X, Y, W, compute.oob.predictions = FALSE, num.trees = 10)
  
  # Join the forests together. 
  big_rf = join_forests(list(c.forest1, c.forest2))
  
  # Result is also a regression_forest of the same class
  expect_true(big_rf$num.trees == (c.forest1$num.trees + c.forest2$num.trees))
  expect_true(is(big_rf, "grf"))
  expect_equal(class(c.forest1), class(big_rf))
})


