#install.packages("~/git/split-relabel/ranger-r-package/ranger", type="source", repos=NULL)

rm(list = ls())

library(ranger)
p = 20
n = 5000

ticks = 1001
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = xvals > 0

X = matrix(2 * runif(n * p) - 1, n, p)
A = rnorm(n)
I = rnorm(n)
W = A + I
Y = 2 * (X[,1] <= 0) * A + (X[,1] > 0) * W + (1 + (sqrt(3) - 1) * (X[,1] > 0)) * rnorm(n)

forest.causal = ranger(Y ~ ., data.frame(X=X, Y=Y, W=W, QUASI_I=I), instrumental = TRUE, num.trees = 100, min.node.size = 200, mtry = p, write.forest = TRUE, status.variable.name="W", instrument.variable.name="QUASI_I", replace=FALSE, sample.fraction=0.2)
preds.causal = predict(forest.causal, data.frame(X=X.test, W=0, QUASI_I=0))$predictions

forest.iv = ranger(Y ~ ., data.frame(X=X, Y=Y, W=W, I=I), instrumental = TRUE, num.trees = 1000, min.node.size = 200, mtry = p, write.forest = TRUE, status.variable.name="W", instrument.variable.name="I", replace=FALSE, sample.fraction=0.2)
preds.iv = predict(forest.iv, data.frame(X=X.test, W=-1, I=-1))$predictions

pdf(paste0("~/git/split-relabel/experiments/IV_forest_causal_splitting.pdf"))
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(xvals), ylim=range(truth,preds.causal,preds.iv, 1.5), xlab="X", ylab="tau")
lines(xvals, truth, lwd = 2, col = 1)
lines(xvals, preds.causal, lwd = 2, col = 4)
lines(xvals, preds.iv, lwd = 2, col = 2)
legend("topleft", c("True Treatment Effect", "IV Forest with Causal Split", "IV Forest with IV Split"), lty=c(1, 1, 1), col=c(1, 4, 2), lwd=2, cex=1.5)
par = pardef
dev.off()