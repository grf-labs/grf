
library(grf)
for(seed in 1:10) {
  set.seed(seed)
  n <- 5000
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-2*X[, 1]))
  W <- rbinom(n, 1, e)
  tau <- 2/e
  Y <- W * tau + rnorm(n)

  e.cc <- 1 / (1+exp(-2*X[,1]))
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- (1-e.cc) / e.cc

  num.trees <- 500
  mse <- function(f) {
      tau.hat <- rep(NA, n)
      tau.hat[cc] <- predict(f)$predictions
      tau.hat[!cc] <- predict(f, X[!cc, ])$predictions
      mean((tau.hat[!cc] - tau[!cc])^2)
  }

  forest.1 <- causal_forest(X[cc, ], Y[cc], W[cc], num.trees = num.trees, seed=seed)
  forest.2 <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights = rep(-1, sum(cc)), num.trees = num.trees, seed=seed)
  forest.3 <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights = -sample.weights[cc], num.trees = num.trees, seed=seed)
  forest.4 <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights  = sample.weights[cc],  num.trees = num.trees, seed=seed)
  print(sapply(list(forest.1, forest.2, forest.3, forest.4), mse))
}

# no weights small-q, no weights, weights but unweighted splits, weights with weighted splits, 
### 14358.09 14356.02 14256.88 14271.42
### 38441.34 38443.80 38149.01 38071.82
### 7090.316 7088.684 7024.002 7027.874
### 6265.315 6271.055 6062.643 6063.425
### 13210.93 13211.71 12821.53 12748.69
### 578373.9 578369.6 578154.7 578186.0
### 7371.337 7374.717 7322.936 7342.382
### 52070.38 52070.28 51926.99 51926.32
### 8799.656 8796.732 8673.229 8675.466
### 10889.09 10890.05 10791.96 10799.49
