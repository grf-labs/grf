library(grf)

num.trees = 500

sample.data = function(seed) {
  set.seed(seed)
  n <- 1000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-X[, 1] + 2 * X[, 2]))
  W <- rbinom(n, 1, e)
  tau <- 2 * (X[, 1] > 0 & X[, 5] > 0) -
      0.5 * (X[, 2] > 0) - 0.5 * (X[, 3] > 0) - 0.5 * (X[, 4] > 0)
  Y <- W * tau + rnorm(n)

  e.cc <- 1 - 0.9 * (X[, 1] > 0 & X[, 5] > 0)
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc
  list(X=X,W=W,Y=Y,tau=tau, e=e, cc=cc, sample.weights=sample.weights)
}

out = sapply(1:50, function(seed) { 
  data = sample.data(seed)
  X = data$X; W=data$W; Y=data$Y; tau=data$tau; e = data$e; cc=data$cc; sample.weights=data$sample.weights;

  mse <- function(f) {
      tau.hat <- rep(NA, length(Y))
      tau.hat[cc] <- predict(f)$predictions
      tau.hat[!cc] <- predict(f, X[!cc, ])$predictions
      mean((tau.hat - tau)^2)
  }

  forest.1 <- causal_forest(X[cc, ], Y[cc], W[cc], num.trees = num.trees, seed=seed)
  forest.2 <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights = -sample.weights[cc], num.trees = num.trees, seed=seed)
  forest.3 <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights = sample.weights[cc], num.trees = num.trees, seed=seed)
  forest.4 <- causal_forest(X[cc, ], Y[cc], W[cc], sample.weights = sample.weights[cc]/sum(sample.weights[cc]), num.trees = num.trees, seed=seed)
  c(mse(forest.2) / mse(forest.1), mse(forest.3) / mse(forest.2), mse(forest.4) / mse(forest.3))
})

rowMeans(out)
# 0.7521323 0.9443328 0.9990205
# 25% improvement using sample weighting in leaves over no weighting
# 5% improvement relative to that using sample weighted splitting
# no meaningful variation to do rescaling of the sample weights
