set.seed(1234)

rm(list = ls())

setwd("~/git/grf/experiments/instrumental_examples")

library(grf)
p = 20
n = 10000

tau = function(x) (2 * as.numeric(x > -1/3))

ticks = 1001
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = tau(xvals)
decoy = truth + 2 * (xvals > 1/3)

X = matrix(2 * runif(n * p) - 1, n, p)
Z = rbinom(n, 1, 1/(1 + exp(X[,2])))
A = rbinom(n, 1, 1/(1 + exp(X[,3])))
W = Z * A
Y = (W - 1/2) * tau(X[,1]) + 3/2 * (2*A - 1) * (X[,1] > 1/3) + 2 * rnorm(n)

forest.causal = causal_forest(X, Y, W, min.node.size = 10, mtry = p)
preds.causal = predict(forest.causal, X.test)$predictions

forest.iv = instrumental_forest(X, Y, W, Z, min.node.size = 10, mtry = p)
preds.iv = predict(forest.iv, X.test)$predictions

pdf("IV_plot_10k_20.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(xvals), ylim=range(truth,preds.causal,preds.iv), xlab="X", ylab="tau")
lines(xvals, truth, lwd = 2, col = 1)
lines(xvals, decoy, lwd = 2, col = 1, lty = 2)
lines(xvals, preds.causal, lwd = 2, col = 4)
lines(xvals, preds.iv, lwd = 2, col = 2)
legend("topleft", c("True Treat. Effect", "Raw Correlation", "Causal Forest", "IV Forest"), lty=c(1, 2, 1, 1), col=c(1, 1, 4, 2), lwd=2, cex=1.5)
par = pardef
dev.off()
