# Check that the following code paths passes under valgrind.
# Usage:
# R -d "valgrind --tool=memcheck --leak-check=full" --vanilla  < test_grf_valgrind.R
library(grf)
p <- 6
n <- 750

X <- round(matrix(runif(n * p), n, p), 2)
W <- rbinom(n, 1, 0.5)
W.factor <- as.factor(W)
D <- rep(1, n)
Y <- (X[, 1] > 0) * (2 * W - 1) + 2 * rnorm(n)
Y.surv <- abs(round(Y, 1))

# Also checks regression_forest
forest.causal <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250, compute.oob.predictions = FALSE)
pp.forest.causal <- predict(forest.causal, X, estimate.variance = TRUE)
pp.forest.causal.oob <- predict(forest.causal, estimate.variance = TRUE)

# Also checks multi_regression_forest and probability_forest
forest.mcausal <- multi_arm_causal_forest(X, Y, W.factor, num.trees = 250, compute.oob.predictions = FALSE)
pp.forest.causal <- predict(forest.mcausal, X, estimate.variance = TRUE)
pp.forest.causal.oob <- predict(forest.mcausal, estimate.variance = TRUE)

# Also checks survival_forest
forest.scausal <- causal_survival_forest(X, Y.surv, W, D, horizon = max(Y.surv), num.trees = 250, compute.oob.predictions = FALSE)
pp.scausal <- predict(forest.scausal, X, estimate.variance = TRUE)
pp.scausal.oob <- predict(forest.scausal, estimate.variance = TRUE)

forest.quantile <- quantile_forest(X, Y, num.trees = 250, compute.oob.predictions = FALSE)
pp.forest.quantile <- predict(forest.quantile, X)
pp.forest.quantile.oob <- predict(forest.quantile)
