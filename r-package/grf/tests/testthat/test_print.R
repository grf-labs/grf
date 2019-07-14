test_that("printing doesn't cause error", {
	p = 4
	n = 50
	i = 2
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,i] > 0))
	D = data.frame(X=X, Y=Y)
	qrf.grad = quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), num.trees = 50)
    capture.output(print(qrf.grad))
})
