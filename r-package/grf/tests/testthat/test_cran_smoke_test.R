library(grf)
set.seed(1234)

test_that("regression forest split frequencies are reasonable", {
	n = 100
	p = 6
	X = matrix(rnorm(n*p), n, p)
	Y = 1000 * (X[,1]) + rnorm(n)
	rrr = regression_forest(X, Y, mtry = p)
	freq = split_frequencies(rrr, 4)
	expect_true(freq[1,1] / sum(freq[1,]) > 1/2)
})

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

    forest.causal = causal_forest(X, Y, W, num.trees = 2000,
                                  ci.group.size = 4, W.hat = 0.5,
                                  compute.oob.predictions = FALSE)
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
