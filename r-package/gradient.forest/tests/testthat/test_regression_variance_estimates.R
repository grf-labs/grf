library(gradient.forest)

set.seed(1234)

test_that("regression variance estimates are positive", {
	p = 2
	n = 400

	ticks = 101
	X.test = matrix(0, ticks, p)
	xvals = seq(-1, 1, length.out = ticks)
	X.test[,1] = xvals
	truth = xvals > 0

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = (X[,1] > 0) + rnorm(n)

	forest = regression.forest(X, Y, sample.fraction = 0.5, num.trees = 10000, ci.group.size = 4)
	preds = predict(forest, X.test, estimate.variance=TRUE)

	expect_true(all(preds$variance.estimate > 0))
})