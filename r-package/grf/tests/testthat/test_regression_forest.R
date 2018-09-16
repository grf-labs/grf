library(grf)

set.seed(1234)

extract_samples <- function(tree) {
  
  # Keep only leaf nodes
  leaf_nodes <- Filter(f = function(x) x$is_leaf, tree$nodes)
  
  # Leaf nodes' 'samples' are estimation samples
  estimation_sample <- unlist(Map(f=function(x) x$samples, leaf_nodes))
  
  # Split = Drawn - Samples
  split_sample <- base::setdiff(tree$drawn_samples, estimation_sample)
  
  return(list(estimation_sample=estimation_sample,
              split_sample=split_sample))
}

test_that("changing honest.fraction behaves as expected", {
  sample_fraction_1 = 0.5
  honesty_fraction_1 = 0.25
  
  sample_fraction_2 = 0.25
  honesty_fraction_2 = 0.1
  
  sample_fraction_3 = 0.25
  honesty_fraction_3 = 0.9
  
  n <- 16
  k <- 10
  X <- matrix(runif(n*k), nrow=n, ncol=k)
  Y <- matrix(runif(n), nrow=n, ncol=1)
  forest_1 <- grf::regression_forest(X, Y, sample.fraction = sample_fraction_1, 
                                   honesty = TRUE, honesty.fraction = honesty_fraction_1)
  samples <- extract_samples(get_tree(forest_1, 1))
  
  expect_equal(length(samples$split_sample), n * sample_fraction_1 * honesty_fraction_1)
  expect_equal(length(samples$estimation_sample), n * sample_fraction_1 * (1 - honesty_fraction_1))
  expect_error(grf::regression_forest(X, Y, sample.fraction = sample_fraction_2, 
                                      honesty = TRUE, honesty.fraction = honesty_fraction_2),
               "The honesty fraction is too close to 1 or 0, as no observations will be sampled.")
  expect_error(grf::regression_forest(X, Y, sample.fraction = sample_fraction_3, 
                                      honesty = TRUE, honesty.fraction = honesty_fraction_3),
               "The honesty fraction is too close to 1 or 0, as no observations will be sampled.")
})

test_that("regression variance estimates are positive", {
	p = 6
	n = 1000

	ticks = 101
	X.test = matrix(0, ticks, p)
	xvals = seq(-1, 1, length.out = ticks)
	X.test[,1] = xvals
	truth = xvals > 0

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = (X[,1] > 0) + 2 * rnorm(n)

	forest = regression_forest(X, Y, num.trees = 1000, ci.group.size = 4)
	preds.oob = predict(forest, estimate.variance=TRUE)
	preds = predict(forest, X.test, estimate.variance=TRUE)
	
	expect_true(all(preds$variance.estimate > 0))
	expect_true(all(preds.oob$variance.estimate > 0))
	
	error = preds$predictions - truth
	expect_true(mean(error^2) < 0.2)
	
	truth.oob = (X[,1] > 0)
	error.oob = preds.oob$predictions - truth.oob
	expect_true(mean(error.oob^2) < 0.2)
	
	Z.oob = error.oob / sqrt(preds.oob$variance.estimate)
	expect_true(mean(abs(Z.oob) > 1) < 0.5)
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

test_that("using a sparse data representation produces the same predictions", {
	dim = 20
	X = diag(rnorm(dim), dim)
	sparse.X = as(X, "dgCMatrix")
	Y = 1000 * (X[,1]) + rnorm(dim)

	forest = regression_forest(X, Y, mtry = dim, seed=10)
	preds = predict(forest, estimate.variance=TRUE)

	sparse.forest = regression_forest(sparse.X, Y, mtry = dim, seed=10)
	sparse.preds = predict(sparse.forest, estimate.variance=TRUE)

	expect_equal(preds$predictions, sparse.preds$predictions)
})

test_that("OOB predictions contain debiased error estimates", {
	p = 6
	n = 10

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = (X[,1] > 0) + 2 * rnorm(n)

	forest = regression_forest(X, Y, num.trees = 1000, ci.group.size = 4)
	preds.oob = predict(forest)

	expect_equal(n, length(preds.oob$debiased.error))
})

test_that("regression forest tuning decreases prediction error", {
	n = 5000
	p = 2

	X = matrix(2 * runif(n * p) - 1, n, p)
	Y = (X[,1] > 0) + rnorm(n)
	X.test = matrix(2 * runif(n * p) - 1, n, p)
	truth = (X.test[,1] > 0)

	forest = regression_forest(X, Y, num.trees = 400, tune.parameters = FALSE)
	preds = predict(forest, X.test)
	error = mean((preds$predictions - truth)^2)
	
	tuned.forest = regression_forest(X, Y, num.trees = 400, tune.parameters= TRUE)
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

  tune.output = tune_regression_forest(X, Y, min.node.size = min.node.size, imbalance.penalty = imbalance.penalty)
  tunable.params = tune.output$params

  expect_equal(as.numeric(tunable.params["min.node.size"]), min.node.size)
  expect_equal(as.numeric(tunable.params["imbalance.penalty"]), imbalance.penalty)
})