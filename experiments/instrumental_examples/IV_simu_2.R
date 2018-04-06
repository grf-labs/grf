set.seed(1234)

rm(list = ls())

setwd("~/git/grf/experiments/instrumental_examples")

library(grf)
p = 20
n = 10000

ticks = 1001
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = xvals > 0

X = matrix(2 * runif(n * p) - 1, n, p)
A = rnorm(n)
Z = rnorm(n)
W = A + Z
Y = 2 * (X[,1] <= 0) * A + (X[,1] > 0) * W + (1 + (sqrt(3) - 1) * (X[,1] > 0)) * rnorm(n)

forest.causal.split = instrumental_forest(X, Y, W, Z, min.node.size = 10, reduced.form.weight = 1, mtry = p)
preds.causal.split = predict(forest.causal.split, X.test)$predictions

forest.iv = instrumental_forest(X, Y, W, Z, min.node.size = 10, reduced.form.weight = 0, mtry = p)
preds.iv = predict(forest.iv, X.test)$predictions

pdf("IV_forest_causal_splitting.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(xvals), ylim=range(truth,preds.causal.split,preds.iv, 1.5), xlab="X", ylab="tau")
lines(xvals, truth, lwd = 2, col = 1)
lines(xvals, preds.causal.split, lwd = 2, col = 4)
lines(xvals, preds.iv, lwd = 2, col = 2)
legend("topleft", c("True Treatment Effect", "IV Forest with Causal Split", "IV Forest with IV Split"), lty=c(1, 1, 1), col=c(1, 4, 2), lwd=2, cex=1.5)
par = pardef
dev.off()
