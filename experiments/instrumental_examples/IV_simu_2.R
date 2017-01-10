rm(list = ls())

setwd("~/git/split-relabel/experiments/instrumental_examples")

library(gradient.forest)
p = 20
n = 5000

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

# TODO: Need to implement Rcpp bindings for causal splitting with IV inference
# Thought... can we add a flag to forest_instrumental that uses causal splitting?
# I can imagine situations where we think W is "almost" unconfounded, and so we use
# causal splitting for more power, but then use IV prediction to be safe.
forest.causal = instrumental.forest(X, Y, W, W, min.node.size = 200)
preds.causal = predict(forest.causal, X.test)

forest.iv = instrumental.forest(X, Y, W, Z, min.node.size = 200)
preds.iv = predict(forest.iv, X.test)

pdf("IV_forest_causal_splitting.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(xvals), ylim=range(truth,preds.causal,preds.iv, 1.5), xlab="X", ylab="tau")
lines(xvals, truth, lwd = 2, col = 1)
lines(xvals, preds.causal, lwd = 2, col = 4)
lines(xvals, preds.iv, lwd = 2, col = 2)
legend("topleft", c("True Treatment Effect", "IV Forest with Causal Split", "IV Forest with IV Split"), lty=c(1, 1, 1), col=c(1, 4, 2), lwd=2, cex=1.5)
par = pardef
dev.off()