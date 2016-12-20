rm(list = ls())

library(quantregForest)
library(ranger)

p = 40
n = 2000

ticks = 1001
X.test = matrix(0, ticks, p)
X.test[,1] = seq(-1, 1, length.out = ticks)
X.test.df = data.frame(X=X.test)

X = matrix(2 * runif(n * p) - 1, n, p)
Y = rnorm(n) * (1 + (X[,1] > 0))
D = data.frame(X=X, Y=Y)

qrf.meinshausen = quantregForest(X, Y, mtry=p, nodesize=10)
preds.meinshausen = predict(qrf.meinshausen, X.test, quantiles = c(0.1, 0.5, 0.9))

qrf.splitrelabel = ranger(Y ~ ., D, quantile = TRUE, num.trees = 1000, min.node.size = 10, mtry = p, write.forest = TRUE, quantiles = c(0.1, 0.5, 0.9))
preds.splitrelabel = predict(qrf.splitrelabel, X.test.df)$predictions

preds.truth = cbind(-qnorm(0.9) * (1 + (X.test[,1] > 0)),
                    0,
                    qnorm(0.9) * (1 + (X.test[,1] > 0))) 

pdf(paste0("~/git/split-relabel/experiments/quantile_plot_spread_n2k_p40.pdf"))

pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim = c(-1, 1), ylim = range(preds.meinshausen, preds.splitrelabel, preds.truth), xlab = "X", ylab = "Y")

lines(X.test[,1], preds.meinshausen[,1], col = 4, lwd = 2, lty = 2)
lines(X.test[,1], preds.meinshausen[,2], col = 4, lwd = 2, lty = 1)
lines(X.test[,1], preds.meinshausen[,3], col = 4, lwd = 2, lty = 2)

lines(X.test[,1], preds.splitrelabel[,1], col = 2, lwd = 2, lty = 2)
lines(X.test[,1], preds.splitrelabel[,2], col = 2, lwd = 2, lty = 1)
lines(X.test[,1], preds.splitrelabel[,3], col = 2, lwd = 2, lty = 2)

lines(X.test[,1], preds.truth[,1], col = 1, lwd = 2, lty = 2)
lines(X.test[,1], preds.truth[,2], col = 1, lwd = 2, lty = 1)
lines(X.test[,1], preds.truth[,3], col = 1, lwd = 2, lty = 2)

legend("topleft", c("truth", "quantregForest", "split-relabel"), lwd = 2, col = c(1, 4, 2), cex=1.5)

par=pardef

dev.off()