set.seed(1234)

rm(list = ls())

setwd("~/git/grf/experiments/quantile_examples")

library(quantregForest)
library(grf)

p = 40
n = 2000

ticks = 1001
X.test = matrix(0, ticks, p)
X.test[,1] = seq(-1, 1, length.out = ticks)
X.test.df = data.frame(X=X.test)

X = matrix(2 * runif(n * p) - 1, n, p)
Y = rnorm(n) * (1 + (X[,1] > 0))
D = data.frame(X=X, Y=Y)

qrf.grad = quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry=p, min.node.size = 10, sample.fraction=0.5)
preds.grad = predict(qrf.grad, X.test.df, quantiles = c(0.1, 0.5, 0.9))

qrf.meinshausen = quantregForest(X, Y, mtry=p, nodesize=10, replace = FALSE, sampsize=ceiling(0.25*n))
preds.meinshausen = predict(qrf.meinshausen, X.test, quantiles = c(0.1, 0.5, 0.9))

preds.truth = cbind(-qnorm(0.9) * (1 + (X.test[,1] > 0)),
                    0,
                    qnorm(0.9) * (1 + (X.test[,1] > 0))) 

pdf("quantile_plot_spread_n2k_p40.pdf")

pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim = c(-1, 1), ylim = range(preds.meinshausen, preds.grad, preds.truth), xlab = "X", ylab = "Y")

lines(X.test[,1], preds.meinshausen[,1], col = 4, lwd = 2, lty = 2)
lines(X.test[,1], preds.meinshausen[,2], col = 4, lwd = 2, lty = 1)
lines(X.test[,1], preds.meinshausen[,3], col = 4, lwd = 2, lty = 2)

lines(X.test[,1], preds.grad[,1], col = 2, lwd = 2, lty = 2)
lines(X.test[,1], preds.grad[,2], col = 2, lwd = 2, lty = 1)
lines(X.test[,1], preds.grad[,3], col = 2, lwd = 2, lty = 2)

lines(X.test[,1], preds.truth[,1], col = 1, lwd = 2, lty = 2)
lines(X.test[,1], preds.truth[,2], col = 1, lwd = 2, lty = 1)
lines(X.test[,1], preds.truth[,3], col = 1, lwd = 2, lty = 2)

legend("topleft", c("truth", "quantregForest", "GRF"), lwd = 2, col = c(1, 4, 2), cex=1.5)

par=pardef

dev.off()
