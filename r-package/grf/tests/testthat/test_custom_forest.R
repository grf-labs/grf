library(grf)

test_that("custom forests behave as expected", {
	p = 40
	n = 500

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,1] > 0))
	D = data.frame(X=X, Y=Y)

	forest = custom_forest(X, Y)
	predictions = predict(forest, X)
	expect_equal(0.0, sum(predictions))
})
