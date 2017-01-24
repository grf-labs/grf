# gradient-forest

![Alt text](documentation/under_construction.jpg?raw=true)

This repository is under active development. To build the R package, go to the `r-package` directory and run the script `build_package.R`. We currently only support Mac OS X (clang 3.3 or higher); linux support is under development.

The software can then be called as follows:
```R
library(gradient.forest)

# generate data
n = 2000; p = 10
X = matrix(rnorm(n*p), n, p)
X.test = matrix(0, 101, p)
X.test[,1] = seq(-2, 2, length.out = 101)

# quantile regression
Y = X[,1] * rnorm(n)
q.forest = quantile.forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
q.hat = predict(q.forest, X.test)
plot(X.test[,1], q.hat[,1], ylim = range(q.hat), xlab = "x", ylab = "y", type = "l")
lines(X.test[,1], q.hat[,2]); lines(X.test[,1], q.hat[,3])
abline(0, qnorm(0.9), lty = 2, col = 2); abline(0, -qnorm(0.9), lty = 2, col = 2)

# treatment effect estimation
W = rbinom(n, 1, 0.5)
Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
tau.forest = causal.forest(X, Y, W)
tau.hat = predict(tau.forest, X.test)
plot(X.test[,1], tau.hat$predictions, ylim = range(tau.hat$predictions, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[,1], pmax(0, X.test[,1]), col = 2, lty = 2)

# add confidence intervals. more trees are now recommended.
tau.forest = causal.forest(X, Y, W, num.trees = 4000)
tau.hat = predict(tau.forest, X.test, estimate.variance = TRUE)
sigma.hat = sqrt(tau.hat$variance.estimates)
plot(X.test[,1], tau.hat$predictions, ylim = range(tau.hat$predictions + 1.96 * sigma.hat, tau.hat$predictions - 1.96 * sigma.hat, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[,1], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[,1], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[,1], pmax(0, X.test[,1]), col = 2, lty = 1)
```
More examples on how to run the functions, including with instrumental variables, can be found in the `experiments` directory.
