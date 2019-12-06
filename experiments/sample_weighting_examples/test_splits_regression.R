library(grf)

out = sapply(1:100, function(seed) { 
    set.seed(seed)
    digits = 1
    n <- 5000
    p <- 2
    X <- round(matrix(rnorm(n * p), n, p),digits)
    Y <- abs(X[, 1]) + 0.1 * rnorm(n)
    e <- 1 / (1 + exp(-2 * X[, 1]))
    sample.weights <- round(1 / e, digits)
    num.trees <- 50
    
    forest.1 <- regression_forest(X, Y, num.trees = num.trees, seed=seed)
    forest.2 <- regression_forest(X, Y, sample.weights = -sample.weights, num.trees=num.trees, seed=seed)
    forest.3 <- regression_forest(X, Y, sample.weights = sample.weights, num.trees=num.trees, seed=seed)
    forest.4 <- regression_forest(X, Y, sample.weights = sample.weights/sum(sample.weights), num.trees=num.trees, seed=seed)
    weighted.mse = function(forest) { mean(sample.weights * (forest$predictions - Y)^2) }
    c(weighted.mse(forest.2)/weighted.mse(forest.1), weighted.mse(forest.3)/weighted.mse(forest.2), weighted.mse(forest.4)/weighted.mse(forest.3)) 
})

rowMeans(out)
# 0.8452165 0.8906758 1.0003165
# 15% improvement using weights within leaves only,
# 10% improvement over that using weights when splitting,
# and no meaningful variation due to rescaling of the sample weights.
# We check the last because when we use weighted spltting, 
# we don't get exactly the same results in the sense that
#    max(abs(predict(forest.3)$predictions - predict(forest.4)$predictions))
#  will be essentially zero. This seems to be rounding error interacting with 
#  the discontinuity of the estimator, not a bug.

