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

test_that("causal forests can split on the last parameter", {
    n = 1000
    p = 6
    X = matrix(rnorm(n*p), n, p)
    W = rbinom(n, 1, 0.5)
    Y = W * (X[,1] + X[,6]) + rnorm(n)
    
    forest = causal_forest(X, Y, W, compute.oob.predictions = FALSE)
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
    ccc = causal_forest(X, Y, W, mtry = p, imbalance.penalty=1.0, stabilize.splits=TRUE, min.node.size=2)
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
    ccc = causal_forest(X, Y, W, mtry = p, imbalance.penalty=1.0, stabilize.splits=FALSE, min.node.size=2)
    split.freq = split_frequencies(ccc, 4)
    expect_true(split.freq[1,p] / sum(split.freq[1,]) > 2/3)
})

test_that("causal forests with a positive imbalance.penalty have reasonable tree depths", {
    n <- 200
    p <- 5
    X <- matrix(rnorm(n * p), n, p)
    W <- rbinom(p=0.5, size=1, n=n)
    Y <- 0.5 * X[,1] * (2*W - 1) + 0.1 * rnorm(n)

    forest = causal_forest(X, Y, W, imbalance.penalty=0.001)
    split.freq = split_frequencies(forest)
    expect_true(sum(split.freq[3,]) > 0)
})

test_that("causal forests with a very small imbalance.penalty behave similarly to unpenalized forests.", {
    n <- 200
    p <- 5
    X <- matrix(rnorm(n * p), n, p)
    W <- rbinom(n, 1, 0.5)
    Y <- X[,1] * (2 * W - 1) + 0.1 * rnorm(n)

    forest = causal_forest(X, Y, W, imbalance.penalty=0.0)
    forest.large.penalty = causal_forest(X, Y, W, imbalance.penalty=100.0)
    forest.small.penalty = causal_forest(X, Y, W, imbalance.penalty=1e-6)

    diff.large.penalty = abs(forest.large.penalty$debiased.error - forest$debiased.error)
    diff.small.penalty = abs(forest.small.penalty$debiased.error - forest$debiased.error)
    expect_true(mean(diff.small.penalty) < 0.10 * mean(diff.large.penalty))
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

test_that("causal forests behave reasonably with small sample size", {
    p = 5
    n = 50
    X = matrix(rnorm(n * p), n, p)
    W = rbinom(n, 1, 0.5)
    tau = 100 * (X[,1] > 0)
    Y = tau * (W - 0.5) + rnorm(n)
    forest = causal_forest(X, Y, W, stabilize.splits = TRUE,
                           min.node.size = 1, mtry = p,
                           compute.oob.predictions = FALSE)
    
    X.test = matrix(rnorm(n * p), n, p)
    tau.test = 100 * (X.test[,1] > 0)
    tau.hat = predict(forest, X.test)$predictions
    expect_true(sqrt(mean((tau.hat - tau.test)^2)) / 100 < 1/3)
})
