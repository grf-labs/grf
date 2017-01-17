library(gradient.forest)

test_that("quantile forests give reasonable results", {
	p = 40
	n = 500

	ticks = 1001
	X.test = matrix(0, ticks, p)
	X.test[,1] = seq(-1, 1, length.out = ticks)
	X.test.df = data.frame(X=X.test)

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,1] > 0))
	D = data.frame(X=X, Y=Y)

	qrf.grad = quantile.forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.632, ci.group.size=1)
	preds.grad = predict(qrf.grad, X.test.df, quantiles = c(0.1, 0.5, 0.9))

	preds.truth = cbind(-qnorm(0.9) * (1 + (X.test[,1] > 0)),
	                    0,
	                    qnorm(0.9) * (1 + (X.test[,1] > 0)))

	mean_error = mean(abs(preds.truth - preds.grad))
	expect_lt(mean_error, 0.5)
})
