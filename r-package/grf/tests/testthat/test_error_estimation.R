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
	
	c(err.debiased.5, err.debiased.20, err.debiased.200)
	c(mse.5, mse.20, mse.200)
	
	expect_equal(err.debiased.200, mse.200, tolerance = 0.01 * sigma)
	
	expect_equal(err.debiased.5, err.debiased.200, tolerance = 0.08 * sigma^2)
	expect_equal(err.debiased.20, err.debiased.200, tolerance = 0.02 * sigma^2)
	
	expect_true(mse.5 - mse.200 >= sigma^2 / 5 / 1.5)
	expect_true(mse.20 - mse.200 >= sigma^2 / 20 / 1.2)
})
