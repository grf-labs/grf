rm(list = ls())

library(quantregForest)

p = 20
n = 1000

ticks = 101
X.test = matrix(0, ticks, p)
X.test[,1] = seq(-1, 1, length.out = ticks)

X = matrix(2 * runif(n * p) - 1, n, p)
Y = rnorm(n) * (1 + (X[,1] > 0))
D = data.frame(X=X, Y=Y)

qrf.meinshausen = quantregForest(X, Y)
preds.meinshausen = predict(qrf.meinshausen, X.test, quantiles = c(0.1, 0.5, 0.9))

qrf.splitrelabel = quantregForest(X, Y)
preds.splitrelabel = predict(qrf.splitrelabel, X.test, quantiles = c(0.1, 0.5, 0.9))

preds.truth = cbind(-qnorm(0.9) * (1 + (X.test[,1] > 0)),
                    0,
                    qnorm(0.9) * (1 + (X.test[,1] > 0))) 

plot(NA, NA, xlim = c(-1, 1), ylim = range(preds.meinshausen, preds.splitrelabel, preds.truth))

lines(X.test[,1], preds.meinshausen[,1], col = 2, lwd = 2, lty = 2)
lines(X.test[,1], preds.meinshausen[,2], col = 2, lwd = 2, lty = 1)
lines(X.test[,1], preds.meinshausen[,3], col = 2, lwd = 2, lty = 2)

lines(X.test[,1], preds.splitrelabel[,1], col = 4, lwd = 2, lty = 2)
lines(X.test[,1], preds.splitrelabel[,2], col = 4, lwd = 2, lty = 1)
lines(X.test[,1], preds.splitrelabel[,3], col = 4, lwd = 2, lty = 2)

lines(X.test[,1], preds.truth[,1], col = 1, lwd = 2, lty = 2)
lines(X.test[,1], preds.truth[,2], col = 1, lwd = 2, lty = 1)
lines(X.test[,1], preds.truth[,3], col = 1, lwd = 2, lty = 2)
