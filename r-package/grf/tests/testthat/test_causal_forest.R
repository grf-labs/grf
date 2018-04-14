library(grf)

set.seed(3141)

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

	forest.causal = causal_forest(X, Y, W, num.trees = 2000, ci.group.size = 4, precompute.nuisance = FALSE)
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

test_that("causal forests have reasonable split frequencies", {
  n = 100
  p = 7
  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.2)
  Y = 1000 * (X[,p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase imbalance.penalty to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc = causal_forest(X, Y, W, mtry = p, imbalance.penalty=0.1, stabilize.splits=TRUE, min.node.size=2)
  split.freq = split_frequencies(ccc, 4)
  expect_true(split.freq[1,p] / sum(split.freq[1,]) > 2/3)
})

test_that("causal forests without stable splitting have reasonable split frequencies", {
  n = 100
  p = 7
  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.2)
  Y = 1000 * (X[,p]) * (2 * W - 1) + rnorm(n)

  # Note that we increase imbalance.penalty to ensure the test reliably passes. Once
  # we add variance corrections, this should no longer be necessary.
  ccc = causal_forest(X, Y, W, mtry = p, imbalance.penalty=0.1, stabilize.splits=FALSE, min.node.size=2)
  split.freq = split_frequencies(ccc, 4)
  expect_true(split.freq[1,p] / sum(split.freq[1,]) > 2/3)
})

test_that("causal forests behave reasonably with a low treatment probability", {
	n = 1000
	p = 5

	X = matrix(rnorm(n * p), n, p)
	W = rbinom(n, 1, 0.1)
	tau = 0.1
	Y = X[,1] + X[,2] + tau * W + rnorm(n)

	forest = causal_forest(X, Y, W, stabilize.splits = TRUE)
	tau.hat = predict(forest)$predictions
	expect_true(sqrt(mean((tau.hat - tau)^2)) < 0.20)
})

test_that("causal forest tuning decreases prediction error", {
	p = 4
	n = 1000

	X = matrix(2 * runif(n * p) - 1, n, p)
	W = rbinom(n, 1, 0.5)
	TAU = 0.1 * (X[,1] > 0)
	Y = TAU * (W  - 1/2) + 2 * rnorm(n)

	forest = causal_forest(X, Y, W, num.trees = 400, min.node.size = 1, tune.parameters = FALSE)
	preds = predict(forest)
	error = mean((preds$predictions - TAU)^2)

	tuned.forest = causal_forest(X, Y, W, num.trees = 400, tune.parameters = TRUE)
	tuned.preds = predict(tuned.forest)
	tuned.error = mean((tuned.preds$predictions - TAU)^2)

	expect_true(tuned.error < error * 0.75)
})

test_that("causal forest tuning only cross-validates null parameters", {
	p = 6
	n = 100

	X = matrix(2 * runif(n * p) - 1, n, p)
	W = rbinom(n, 1, 0.5)
	TAU = 2 * (X[,1] > 0)
	Y = TAU * (W  - 1/2) + 2 * rnorm(n)

	min.node.size = 5
	imbalance.penalty = 0.42

  tune.output = tune_causal_forest(X, Y, W, min.node.size = min.node.size, imbalance.penalty = imbalance.penalty)
  tunable.params = tune.output$params

  expect_equal(as.numeric(tunable.params["min.node.size"]), min.node.size)
  expect_equal(as.numeric(tunable.params["imbalance.penalty"]), imbalance.penalty)
})

test_that("causal forests behave reasonably with small sample size", {
  p = 5
  n = 50
  X = matrix(rnorm(n * p), n, p)
  W = rbinom(n, 1, 0.5)
  tau = 100 * (X[,1] > 0)
  Y = tau * (W - 0.5) + rnorm(n)
  forest = causal_forest(X, Y, W, stabilize.splits = TRUE, min.node.size = 1, mtry = p)
  tau.hat = predict(forest)$predictions
  expect_true(sqrt(mean((tau.hat - tau)^2)) / 100 < 1/4)
})
