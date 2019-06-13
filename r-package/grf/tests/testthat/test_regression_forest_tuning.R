library(grf)

seed <- 1000
set.seed(seed)

test_that("regression forest tuning decreases prediction error", {
	n = 5000
	p = 2

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = (X[,1] > 0) + rnorm(n)
	X.test = matrix(2 * runif(n * p) - 1, n, p)
	truth = (X.test[,1] > 0)

	forest = regression_forest(X, Y, num.trees = 400, tune.parameters = FALSE, seed = seed)
	preds = predict(forest, X.test)
	error = mean((preds$predictions - truth)^2)

	tuned.forest = regression_forest(X, Y, num.trees = 400, tune.parameters = TRUE, seed = seed)
	tuned.preds = predict(tuned.forest, X.test)
	tuned.error = mean((tuned.preds$predictions - truth)^2)

	expect_true(tuned.error < error * 0.75)
})

test_that("regression forest tuning only cross-validates null parameters", {
	n = 5000
	p = 2

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = (X[,1] > 0) + rnorm(n)
	X.test = matrix(2 * runif(n * p) - 1, n, p)
	truth = (X.test[,1] > 0)

	min.node.size = 42
	imbalance.penalty = 0.42

    tune.output = tune_regression_forest(X, Y, min.node.size = min.node.size,
                                         imbalance.penalty = imbalance.penalty, seed = seed)
    tunable.params = tune.output$params

    expect_equal(as.numeric(tunable.params["min.node.size"]), min.node.size)
    expect_equal(as.numeric(tunable.params["imbalance.penalty"]), imbalance.penalty)
})
