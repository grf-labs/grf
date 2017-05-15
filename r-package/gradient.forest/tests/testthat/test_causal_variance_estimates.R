library(gradient.forest)

set.seed(1234)

test_that("causal variance estimates are positive", {
	p = 2
	n = 400

	ticks = 101
	X.test = matrix(0, ticks, p)
	xvals = seq(-1, 1, length.out = ticks)
	X.test[,1] = xvals
	truth = xvals > 0

	X = matrix(2 * runif(n * p) - 1, n, p)
	W = rbinom(n, 1, 0.5)
	Y = (X[,1] > 0) * (2 * W  - 1) + rnorm(n)

	forest.causal = causal.forest(X, Y, W, num.trees = 1000, ci.group.size = 4)
	preds.causal = predict(forest.causal, X.test, estimate.variance=TRUE)

	expect_true(all(preds.causal$variance.estimate > 0))
})