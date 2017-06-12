library(gradient.forest)

# TODO(swager): debug why this doesn't always return the right index!
# test_that("computing split frequencies gives reasonable results", {
# 	p = 40
# 	n = 500

# 	i = 5
# 	X = matrix(2 * runif(n * p) - 1, n, p)
# 	Y = rnorm(n) * (1 + (X[,i] > 0))
# 	D = data.frame(X=X, Y=Y)

# 	qrf.grad = quantile.forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.632, ci.group.size=1)
# 	split.frequencies = compute_split_frequencies(qrf.grad, 10)

# 	expect_equal(which.max(split.frequencies[1,]), i)
# })

test_that("examining a tree gives reasonable results", {
	p = 40
	n = 500

	i = 5
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + (X[,i] > 0))
	D = data.frame(X=X, Y=Y)

	qrf.grad = quantile.forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.632, ci.group.size=1)
	quantile.tree = examine_tree(qrf.grad, 500)

	expect_lt(length(quantile.tree$nodes), n)

    split.vars = unlist(sapply(quantile.tree$nodes, function(node)node$split_variable))
	expect_true(all(split.vars >= 0 && split.vars <= 40))
})
