library(grf)

p <- 6
n <- 1000

ticks <- 101
X.test <- matrix(0, ticks, p)
xvals <- seq(-1, 1, length.out = ticks)
X.test[, 1] <- xvals
truth <- 2 * (xvals > 0)

X <- matrix(2 * runif(n * p) - 1, n, p)
W <- rbinom(n, 1, 0.5)
Y <- (X[, 1] > 0) * (2 * W - 1) + 2 * rnorm(n)

forest.causal <- causal_forest(X, Y, W,
  num.trees = 2000,
  ci.group.size = 4, W.hat = 0.5,
  compute.oob.predictions = FALSE
)
preds.causal.oob <- predict(forest.causal, estimate.variance = TRUE)
preds.causal <- predict(forest.causal, X.test, estimate.variance = TRUE)
preds.causal2 <- predict(forest.causal, X, estimate.variance = TRUE)

forest.mcausal <- multi_arm_causal_forest(X, Y, as.factor(W), compute.oob.predictions = FALSE)
preds.mcausal.oob <- predict(forest.mcausal, estimate.variance = TRUE)

forest.scausal <- causal_survival_forest(X, round(Y, 1), W, rep(1, n),
                                         W.hat = 0.5, compute.oob.predictions = FALSE)
preds.scausal.oob <- predict(forest.scausal, estimate.variance = TRUE)
