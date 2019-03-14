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
    
    tune.output = tune_causal_forest(X, Y, W, 0, 0.5,
                                     min.node.size = min.node.size,
                                     imbalance.penalty = imbalance.penalty)
    tunable.params = tune.output$params
    
    expect_equal(as.numeric(tunable.params["min.node.size"]), min.node.size)
    expect_equal(as.numeric(tunable.params["imbalance.penalty"]), imbalance.penalty)
})
