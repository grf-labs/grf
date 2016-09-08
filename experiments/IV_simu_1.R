#install.packages("~/git/split-relabel/ranger-r-package/ranger", type="source", repos=NULL)

rm(list = ls())

library(ranger)
p = 20
n = 10000

tau = function(x) (2 * as.numeric(x > -1/3))

ticks = 1001
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = tau(xvals)
decoy = truth + 2 * (xvals > 1/3)

replicate(5, {
X = matrix(2 * runif(n * p) - 1, n, p)
I = rbinom(n, 1, 1/(1 + exp(X[,2])))
A = rbinom(n, 1, 1/(1 + exp(X[,3])))
W = I * A
Y = (W - 1/2) * tau(X[,1]) + 3/2 * (2*A - 1) * (X[,1] > 1/3) + 2 * rnorm(n)

forest.causal = ranger(Y ~ ., data.frame(X=X, Y=Y, W=W), causal = TRUE, num.trees = 1000, min.node.size = 10, mtry = p, write.forest = TRUE, status.variable.name="W")
preds.causal = predict(forest.causal, data.frame(X=X.test, W=-1))$predictions

forest.iv = ranger(Y ~ ., data.frame(X=X, Y=Y, W=W, I=I), instrumental = TRUE, num.trees = 1000, min.node.size = 200, mtry = p, write.forest = TRUE, status.variable.name="W", instrument.variable.name="I", replace=FALSE, sample.fraction=0.2)
preds.iv = predict(forest.iv, data.frame(X=X.test, W=-1, I=-1))$predictions

rnd = sample.int(1000, 1)

pdf(paste0("~/git/split-relabel/experiments/IV_plot_10k_20_", rnd, ".pdf"))
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(xvals), ylim=range(truth,preds.causal,preds.iv), xlab="X", ylab="tau")
lines(xvals, truth, lwd = 2, col = 1)
lines(xvals, decoy, lwd = 2, col = 1, lty = 2)
lines(xvals, preds.causal, lwd = 2, col = 4)
lines(xvals, preds.iv, lwd = 2, col = 2)
legend("topleft", c("True Treat. Effect", "Raw Correlation", "Causal Forest", "IV Forest"), lty=c(1, 2, 1, 1), col=c(1, 1, 4, 2), lwd=2, cex=1.5)
par = pardef
dev.off()
})