set.seed(2)

rm(list = ls())

setwd("~/git/grf/experiments/instrumental_examples")

library(grf)
p = 20
n = 4000
num.trees = 10000

ticks = 1001
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = xvals > 0

res = lapply(1:4, function(iter){
	
print(iter)
	
X = matrix(2 * runif(n * p) - 1, n, p)
A = rnorm(n)
Z = rnorm(n)
W = A + Z
Y = 2 * (X[,1] <= 0) * A + (X[,1] > 0) * W + (1 + (sqrt(3) - 1) * (X[,1] > 0)) * rnorm(n)


forest.iv = instrumental_forest(X, Y, W, Z, num.trees = num.trees, ci.group.size = 20, mtry = p)
preds.iv = predict(forest.iv, X.test, estimate.variance=TRUE)
tau.hat = preds.iv$predictions
sigma.hat = sqrt(preds.iv$variance.estimate)

cbind(tau.hat, tau.hat - 1.96 * sigma.hat, tau.hat + 1.96 * sigma.hat)
})

rr = range(sapply(res, range))

pdf("IV_plot_with_CI.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(xvals), ylim=range(truth,rr), xlab="X", ylab="tau")
lines(xvals, truth, lwd = 2, col = 2)
for (iter in 1:4) {
lines(xvals, res[[iter]][,1], lwd = 1, col = 1)
lines(xvals, res[[iter]][,2], lwd = 1, col = 1, lty = 2)
lines(xvals, res[[iter]][,3], lwd = 1, col = 1, lty = 2)
}
legend("topleft", c("True Effect", "Point Estimate", "Confidence Int."), lty=c(1, 1, 2), col=c(2, 1, 1), lwd=c(2, 1, 1), cex=1.5)
par = pardef
dev.off()
