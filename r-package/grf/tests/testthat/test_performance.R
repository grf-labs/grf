library(grf)

### REGRESSION TEST ###
# n = 1000
# p = 100

# ticks = 101
# X.test = matrix(0, ticks, p)
# xvals = seq(-1, 1, length.out = ticks)
# X.test[,1] = xvals
# truth = xvals > 0

# X = matrix(2 * runif(n * p) - 1, n, p)
# Y = (X[,1] > 0) + 2 * rnorm(n)

# forest = regression_forest(X, Y, num.trees = 1000, ci.group.size = 4)
# preds.oob = predict(forest, estimate.variance=TRUE)
# preds = predict(forest, X.test, estimate.variance=TRUE)

### CAUSAL TEST ###
p = 6
n = 10000

X = matrix(2 * runif(n * p) - 1, n, p)
W = rbinom(n, 1, 0.5)
Y = (X[,1] > 0) * (2 * W  - 1) + 2 * rnorm(n)

forest.causal = causal_forest(X, Y, W, num.trees = 500, ci.group.size = 4, sample.fraction = 0.10, precompute.nuisance = TRUE)

# Sys.sleep(5)

ticks = 101
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = 2 * (xvals > 0)

preds.causal.oob = predict(forest.causal, estimate.variance=TRUE)
preds.causal = predict(forest.causal, X.test, estimate.variance=TRUE)
