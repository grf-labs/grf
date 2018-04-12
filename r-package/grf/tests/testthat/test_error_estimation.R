library(grf)

set.seed(1234)

test_that("regression error estimates are reasonable", {
	p = 3
	n = 2000
	sigma = 1

	X = matrix(2 * runif(n * p) - 1, n, p)
	MU = 0.1 * (X[,1] > 0)
	Y = MU + sigma * rnorm(n)

	forest.5 = regression_forest(X, Y, num.trees = 5, ci.group.size = 1, min.node.size = 1)
	pred.5 = predict(forest.5)
	err.debiased.5 = mean(pred.5$debiased.error, na.rm = TRUE)
	mse.5 = mean((pred.5$predictions - Y)^2, na.rm = TRUE)
	
	forest.20 = regression_forest(X, Y, num.trees = 20, ci.group.size = 1, min.node.size = 1)
	pred.20 = predict(forest.20)
	err.debiased.20 = mean(pred.20$debiased.error, na.rm = TRUE)
	mse.20 = mean((pred.20$predictions - Y)^2, na.rm = TRUE)
	
	forest.200 = regression_forest(X, Y, num.trees = 200, ci.group.size = 1)
	pred.200 = predict(forest.200)
	err.debiased.200 = mean(pred.200$debiased.error, na.rm = TRUE)
	mse.200 = mean((pred.200$predictions - Y)^2, na.rm = TRUE)
	
	# c(err.debiased.5, err.debiased.20, err.debiased.200)
	# c(mse.5, mse.20, mse.200)
	
	expect_equal(err.debiased.200, mse.200, tolerance = 0.01 * sigma^2)
	expect_equal(err.debiased.5, err.debiased.200, tolerance = 0.08 * sigma^2)
	expect_equal(err.debiased.20, err.debiased.200, tolerance = 0.02 * sigma^2)
	
	expect_true(mse.5 - mse.200 >= sigma^2 / 5 / 1.5)
	expect_true(mse.20 - mse.200 >= sigma^2 / 20 / 1.2)
})

test_that("causal error estimates are reasonable", {
  p = 3
  n = 2000
  sigma = 0.1
  
  X = matrix(2 * runif(n * p) - 1, n, p)
  W = rbinom(n, 1, 0.1)
  TAU = (X[,1] > 0)
  Y = TAU * (W - 1/2) + sigma * rnorm(n)
  
  W.forest = regression_forest(X, W, num.trees = 500, sample.fraction = 0.2)
  W.hat = predict(W.forest)$predictions
  
  Y.forest = regression_forest(X, Y, num.trees = 500, sample.fraction = 0.2)
  Y.hat = predict(Y.forest)$predictions
  
  Y.resid = Y - Y.hat
  W.resid = W - W.hat
  
  cf.10 = causal_forest(X, Y.resid, W.resid, precompute.nuisance = FALSE, num.trees = 10,
                        min.node.size = 1, stabilize.splits = TRUE)
  tau.hat.10 = predict(cf.10)
  raw.10 = mean((Y.resid - tau.hat.10$predictions * W.resid)^2, na.rm = TRUE)
  err.10 = mean(tau.hat.10$debiased.error, na.rm = TRUE)
  
  cf.20 = causal_forest(X, Y.resid, W.resid, precompute.nuisance = FALSE, num.trees = 20,
                        min.node.size = 1, stabilize.splits = TRUE)
  tau.hat.20 = predict(cf.20)
  raw.20 = mean((Y.resid - tau.hat.20$predictions * W.resid)^2, na.rm = TRUE)
  err.20 = mean(tau.hat.20$debiased.error, na.rm = TRUE)
  
  cf.400 = causal_forest(X, Y.resid, W.resid, precompute.nuisance = FALSE, num.trees = 400,
                         min.node.size = 1, stabilize.splits = TRUE)
  tau.hat.400 = predict(cf.400)
  raw.400 = mean((Y.resid - tau.hat.400$predictions * W.resid)^2)
  err.400 =  mean(tau.hat.400$debiased.error, na.rm = TRUE)
  
  # c(err.10, err.20, err.400) / sigma^2
  # c(raw.10, raw.20, raw.400) / sigma^2
  
  expect_equal(err.400, raw.400, tolerance = 0.1 * sigma^2)
  expect_equal(err.10, err.400, tolerance = 0.5 * sigma^2)
  expect_equal(err.20, err.400, tolerance = 0.3 * sigma^2)
})
