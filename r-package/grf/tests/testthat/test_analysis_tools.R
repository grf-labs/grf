library(grf)

set.seed(1234)

test_that("quantile forest split frequencies are reasonable", {
	p = 10
	n = 500
	i = 5
	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = rnorm(n) * (1 + 100 * (X[,i] > 0))
	qrf = quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.632)
	split.frequencies = split_frequencies(qrf, 4)
	expect_true(split.frequencies[1,i]/sum(split.frequencies[1,]) > 1/2)
})

test_that("regression forest split frequencies are reasonable", {
  n = 100
  p = 6
  X = matrix(rnorm(n*p), n, p)
  Y = 1000 * (X[,1]) + rnorm(n)
  rrr = regression_forest(X, Y, mtry = p)
  freq = split_frequencies(rrr, 4)
  expect_true(freq[1,1] / sum(freq[1,]) > 1/2)
})

test_that("causal forest split frequencies are reasonable", {
  n = 100
  p = 7
  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.2)
  Y = 1000 * (X[,p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase alpha to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc = causal_forest(X, Y, W, mtry = p, alpha=0.05)
  freq = split_frequencies(ccc, 4)
  expect_true(freq[1,p] / sum(freq[1,]) > 2/3)
})

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
