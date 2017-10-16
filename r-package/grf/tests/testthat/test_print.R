test_that("basic printing is successful", {
	p = 4
	n = 50
	i = 2
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,i] > 0))
	D = data.frame(X=X, Y=Y)
	q.forest = quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), num.trees = 50)
    capture.output(print(q.forest))
})

test_that("printing a forest with one regressor is successful", {
	n = 50; p = 1
	X = matrix(rnorm(n * p), n, p)
	Y = X * rnorm(n)
	r.forest = regression_forest(X, Y)
	capture.output(print(r.forest))
})
