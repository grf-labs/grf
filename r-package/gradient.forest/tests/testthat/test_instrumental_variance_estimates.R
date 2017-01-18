library(gradient.forest)

set.seed(1234)

test_that("instrumental variance estimates are reasonable", {
	
	p = 20
	n = 1000

	ticks = 101
	X.test = matrix(0, ticks, p)
	xvals = seq(-1, 1, length.out = ticks)
	X.test[,1] = xvals
	truth = xvals > 0

	X = matrix(2 * runif(n * p) - 1, n, p)
	A = rnorm(n)
	Z = rnorm(n)
	W = A + Z
	Y = 2 * (X[,1] <= 0) * A + (X[,1] > 0) * W + (1 + (sqrt(3) - 1) * (X[,1] > 0)) * rnorm(n)

	forest.iv = instrumental.forest(X, Y, W, Z, min.node.size = 5, split.regularization = 0)
	preds.iv = predict(forest.iv, X.test)
})