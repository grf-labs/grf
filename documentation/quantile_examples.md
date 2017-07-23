```R
library(grf)

# Generate data.
n = 2000; p = 20
X = matrix(rnorm(n*p), n, p)
X.test = matrix(0, 101, p)
X.test[,1] = seq(-2, 2, length.out = 101)
Y = X[,1] * rnorm(n)

# Train a quantile forest.
q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))

# Make predictions.
q.hat = predict(q.forest, X.test)
plot(X.test[,1], q.hat[,1], ylim = range(q.hat), xlab = "x", ylab = "y", type = "l")
lines(X.test[,1], q.hat[,2]); lines(X.test[,1], q.hat[,3])
abline(0, qnorm(0.9), lty = 2, col = 2); abline(0, -qnorm(0.9), lty = 2, col = 2)

# Make predictions for different quantiles than those used in training.
q.hat = predict(q.forest, X.test, quantiles=c(0.1, 0.9))
plot(X.test[,1], q.hat[,1], ylim = range(q.hat), xlab = "x", ylab = "y", type = "l")
lines(X.test[,1], q.hat[,2])
abline(0, qnorm(0.9), lty = 2, col = 2); abline(0, -qnorm(0.9), lty = 2, col = 2)

# Train a quantile forest using regression splitting instead of quantile-based splits.
# Note that this emulates the approach in Meinshausen (2006).
q.forest = quantile_forest(X, Y, regression.splitting=TRUE)

# Make predictions for the desired quantiles.
q.hat = predict(q.forest, X.test, quantiles=c(0.1, 0.5, 0.9))
plot(X.test[,1], q.hat[,1], ylim = range(q.hat), xlab = "x", ylab = "y", type = "l")
lines(X.test[,1], q.hat[,2])
abline(0, qnorm(0.9), lty = 2, col = 2); abline(0, -qnorm(0.9), lty = 2, col = 2)
```