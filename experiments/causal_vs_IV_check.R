p = 10
n = 1000


tau = function(x) (2 * as.numeric(x > 0))

ticks = 1001
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = tau(xvals)

X = matrix(2 * runif(n * p) - 1, n, p)
W = rbinom(n, 1, 1/(1 + exp(X[,3])))
Y = (W - 1/2) * tau(X[,1]) + X[,2] + rnorm(n)

forest.causal = ranger(Y ~ ., data.frame(X=X, Y=Y, W=W), causal = TRUE, num.trees = 1000, min.node.size = 10, mtry = p, write.forest = TRUE, status.variable.name="W")
preds.causal = predict(forest.causal, data.frame(X=X.test, W=-1))$predictions

forest.iv = ranger(Y ~ ., data.frame(X=X, Y=Y, W=W, I=W), instrumental = TRUE, num.trees = 1000, min.node.size = 10, mtry = p, write.forest = TRUE, status.variable.name="W", instrument.variable.name="I")
preds.iv = predict(forest.iv, data.frame(X=X.test, W=-1, I=-1))$predictions

plot(NA, NA, xlim=range(xvals), ylim=range(truth,preds.causal,preds.iv), xlab="X", ylab="tau")
lines(xvals, truth, lwd = 2, col = 1)
lines(xvals, preds.causal, lwd = 2, col = 4)
lines(xvals, preds.iv, lwd = 2, col = 2)
legend("topleft", c("True Treat. Effect", "Causal Forest", "IV Forest"), lty=c(1), col=c(1, 4, 2), lwd=2)