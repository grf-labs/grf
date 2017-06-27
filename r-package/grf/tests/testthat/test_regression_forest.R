library(grf)

set.seed(1234)

test_that("regression variance estimates are positive", {
	p = 6
	n = 1000

	ticks = 101
	X.test = matrix(0, ticks, p)
	xvals = seq(-1, 1, length.out = ticks)
	X.test[,1] = xvals
	truth = xvals > 0

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = (X[,1] > 0) + 2 * rnorm(n)

	forest = regression_forest(X, Y, num.trees = 1000, ci.group.size = 4)
	preds.oob = predict(forest, estimate.variance=TRUE)
	preds = predict(forest, X.test, estimate.variance=TRUE)
	
	expect_true(all(preds$variance.estimate > 0))
	expect_true(all(preds.oob$variance.estimate > 0))
	
	error = preds$predictions - truth
	expect_true(mean(error^2) < 0.2)
	
	truth.oob = (X[,1] > 0)
	error.oob = preds.oob$predictions - truth.oob
	expect_true(mean(error.oob^2) < 0.2)
	
	Z.oob = error.oob / sqrt(preds.oob$variance.estimate)
	expect_true(mean(abs(Z.oob) > 1) < 0.5)
})

test_that("regression forest split frequencies are reasonable", {
  n = 100
  p = 6
  X = matrix(rnorm(n*p), n, p)
  Y = 1000 * (X[,1]) + rnorm(n)
  rrr = regression_forest(X, Y, mtry = p)
  freq = split_frequencies(rrr, 4)
  expect_true(freq[1,1] / sum(freq[1,]) > 1/2)
})
