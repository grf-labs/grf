library(gradient.forest)

set.seed(1234)

test_that("causal variance estimates are reasonable", {

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

	forest.causal = instrumental.forest(X, Y, W, W, sample.fraction = 0.5, num.trees = 10000, split.regularization = 0, ci.group.size = 4)
	preds.causal = predict(forest.causal, X.test, estimate.variance=TRUE)

	preds.med = preds.causal[[1]]
	preds.low = preds.causal[[1]] - 1.96 * sqrt(pmax(0, preds.causal[[2]]))
	preds.high = preds.causal[[1]] + 1.96 * sqrt(pmax(0, preds.causal[[2]]))

	plot(X.test[,1], preds.med, ylim = range(preds.med, preds.low, preds.high))
    lines(X.test[,1], preds.low)
	lines(X.test[,1], preds.high)

})