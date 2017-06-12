test_that("computing split frequencies gives reasonable results", {
	p = 40
	n = 500

	i = 5
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,i] > 0))
	D = data.frame(X=X, Y=Y)

	qrf.grad = quantile.forest(X, Y, quantiles = c(0.1, 0.5, 0.9))
	invisible(print(qrf.grad))
})