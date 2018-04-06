library(grf)

set.seed(4321)

test_that("causal forests give reasonable estimates", {
	p = 6
	n = 1000

	ticks = 101
	X.test = matrix(0, ticks, p)
	xvals = seq(-1, 1, length.out = ticks)
	X.test[,1] = xvals
	truth = 2 * (xvals > 0)

	X = matrix(2 * runif(n * p) - 1, n, p)
	W = rbinom(n, 1, 0.5)
	Y = (X[,1] > 0) * (2 * W  - 1) + 2 * rnorm(n)

	forest.causal = causal_forest(X, Y, W, num.trees = 2000, ci.group.size = 4)
	preds.causal.oob = predict(forest.causal, estimate.variance=TRUE)
	preds.causal = predict(forest.causal, X.test, estimate.variance=TRUE)

	expect_true(all(preds.causal$variance.estimate > 0))
	expect_true(all(preds.causal.oob$variance.estimate > 0))
	
	error = preds.causal$predictions - truth
	expect_true(mean(error^2) < 0.5)
	
	truth.oob = 2 * (X[,1] > 0)
	error.oob = preds.causal.oob$predictions - truth.oob
	expect_true(mean(error.oob^2) < 0.5)
	
	Z.oob = error.oob / sqrt(preds.causal.oob$variance.estimate)
	expect_true(mean(abs(Z.oob) > 1) < 0.5)
})

test_that("causal forests can split on the last parameter", {
	n = 1000
	p = 6
	X = matrix(rnorm(n*p), n, p)
	W = rbinom(n, 1, 0.5)
	Y = W * (X[,1] + X[,6]) + rnorm(n)
	 
	forest = causal_forest(X, Y, W)
	split.freq = split_frequencies(forest, 10)

 	expect_gt(sum(split.freq[,6]), 0)
})

test_that("causal forest split frequencies are reasonable", {
  n = 100
  p = 7
  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.2)
  Y = 1000 * (X[,p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase imbalance.penalty to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc = causal_forest(X, Y, W, mtry = p, imbalance.penalty=0.1)
  split.freq = split_frequencies(ccc, 4)
  expect_true(split.freq[1,p] / sum(split.freq[1,]) > 2/3)
})

test_that("causal forests behave reasonably with a low treatment probability", {
	p = 20
	n = 1000
	X = matrix(rnorm(n * p), n, p)
	W = rbinom(n, 1, 0.01)
	tau = 0.1
	Y = X[,1] + X[,2] + tau * W + rnorm(n)

	ff = causal_forest(X, Y, W, stabilize.splits = TRUE, min.node.size = 30)
	tau.hat = predict(ff)$predictions

	print(sqrt(mean((tau.hat - tau)^2)))
	expect_true(sqrt(mean((tau.hat - tau)^2)) < 0.20)
})
