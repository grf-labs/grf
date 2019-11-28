library(grf)

sapply(1:10, function(seed) { 
    set.seed(seed)
    n <- 5000
    p <- 2
    X <- matrix(rnorm(n * p), n, p)
    Y <- abs(X[, 1]) + 0.1 * rnorm(n)
    e <- 1 / (1 + exp(-2 * X[, 1]))
    sample.weights <- 1 / e
    num.trees <- 500

    forest.1 <- regression_forest(X, Y, num.trees = num.trees, seed=seed)
    forest.2 <- regression_forest(X, Y, sample.weights = -sample.weights, num.trees=num.trees, seed=seed)
    forest.3 <- regression_forest(X, Y, sample.weights =  sample.weights, num.trees=num.trees, seed=seed)
    weighted.mse = function(forest) { mean(sample.weights * (forest$predictions - Y)^2) }
    print(sapply(list(forest.1, forest.2, forest.3), weighted.mse)) 
})

# small-q unweighted, weights but unweighted splits, weights with weighted splits
### 0.3624662 0.3142132 0.3050226
### 1.715285  1.297062  1.079989
### 0.3226483 0.2922827 0.2786545
### 0.2606093 0.1946755 0.1828769
### 0.4085726 0.3301638 0.2959058
### 15.66919  14.69359  13.96520
### 0.3192164 0.2726220 0.2167435
### 1.2060341 0.8869051 0.8183573
### 0.3395752 0.2765931 0.2580846
### 0.3355356 0.2905805 0.2407564
