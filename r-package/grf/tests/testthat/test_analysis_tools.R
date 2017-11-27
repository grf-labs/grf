library(grf)

set.seed(1234)

test_that("examining a tree gives reasonable results", {
	p = 40
	n = 500

	i = 5
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,i] > 0))
	D = data.frame(X=X, Y=Y)

	qrf = quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.632)
	quantile.tree = get_tree(qrf, 500)

	expect_lt(length(quantile.tree$nodes), n)

  	split.vars = unlist(sapply(quantile.tree$nodes, function(node)node$split_variable))
	expect_true(all(split.vars >= 0 && split.vars <= 40))
})

test_that("tree indexes are one-indexed", {
	p = 40
	n = 500

	i = 5
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,i] > 0))
	D = data.frame(X=X, Y=Y)

	qrf = quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.632)

	expect_error(get_tree(qrf, 0))
	get_tree(qrf, n)
})

test_that("calculating variable importance gives reasonable results", {
	p = 40
	n = 500

	i = 5
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,i] > 0))
	D = data.frame(X=X, Y=Y)

	qrf = quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.632)

	var.importance = variable_importance(qrf)

	expect_equal(length(var.importance), p)
	expect_equal(which.max(var.importance), i)
})
