rm(list = ls())

setwd("~/git/grf/experiments/instrumental_stress_test")

library(grf)
source("../baselines.R")

#
#
# PROBLEM CONFIGURATION
#
#

alpha.tau = 1
alpha.mu.vals = c(0, 3)
sigma = 1
confounding.vals = c(0, 1)
axis.aligned.vals = c(TRUE, FALSE)

k.tau.vals = c(2, 4)
k.mu = 2
p.vals = c(10, 20)

n.vals = c(1000, 2000)
n.test = 1000

REPS = 100

for (alpha.mu in alpha.mu.vals) {
for (confounding in confounding.vals) {
for (axis.aligned in axis.aligned.vals) {
for (k.tau in k.tau.vals) {
for (p in p.vals) {
for (n in n.vals) {

p = max(p, k.tau + k.mu)

res.all = replicate(REPS, {

#
#
# DATA GENERATION
#
#

X = matrix(rnorm(n * p), n, p)
eps = rnorm(n)
Z = rbinom(n, 1, 2/3)
filter = rbinom(n, 1, 1/(1 + exp(-confounding * eps)))
W = Z * filter

if(axis.aligned) {
  tau = alpha.tau * apply(X[,1:k.tau], 1, function(xx) sum(pmax(0, xx)))
  mu = alpha.mu *  apply(X[,k.tau + 1:k.mu], 1, function(xx) sum(pmax(0, xx)))
} else {
  tau = alpha.tau * pmax(0, rowSums(X[,1:k.tau]))
  mu = alpha.mu * pmax(0, rowSums(X[,k.tau + 1:k.mu]))
}

Y = (2 * W - 1) / 2 * tau + mu + eps

X.test = matrix(rnorm(n.test * p), n.test, p)
if(axis.aligned) {
  tau.true = alpha.tau *  apply(X.test[,1:k.tau], 1, function(xx) sum(pmax(0, xx)))
} else {
  tau.true = alpha.tau * pmax(0, rowSums(X.test[,1:k.tau]))
}

#
#
# RUN EXPERIMENT
#
#

forest.orthog = instrumental_forest(X, Y, W, Z, precompute.nuisance = TRUE)
tau.forest.orthog = predict(forest.orthog, newdata = X.test)$predictions

forest.plain = instrumental_forest(X, Y, W, Z, precompute.nuisance = FALSE)
tau.forest.plain = predict(forest.plain, newdata = X.test)$predictions

tau.iv.noint = iv.series(X, Y, W, Z, X.test, interact = FALSE)
tau.iv.int = iv.series(X, Y, W, Z, X.test, interact = TRUE)

tau.knn.10 = iv.knn(X, Y, W, Z, X.test, k = 10)
tau.knn.30 = iv.knn(X, Y, W, Z, X.test, k = 30)
tau.knn.100 = iv.knn(X, Y, W, Z, X.test, k = 100)
tau.knn.300 = iv.knn(X, Y, W, Z, X.test, k = 300)

mse = function(tt) mean((tt - tau.true)^2)

res = c(forest.orthog=mse(tau.forest.orthog),
        forest.plain=mse(tau.forest.plain),
        iv.noint=mse(tau.iv.noint),
        iv.int=mse(tau.iv.int),
        knn.10=mse(tau.knn.10),
        knn.30=mse(tau.knn.30),
        knn.100=mse(tau.knn.100),
        knn.300=mse(tau.knn.300))

print(round(res, 2))

res

})

print(paste("CONFIG", alpha.tau, alpha.mu, sigma, confounding, axis.aligned, k.tau, k.mu, p, n, sep="-"))
print(rowMeans(res.all))
print(sqrt(apply(res.all, 1, var)/REPS))
print(apply(res.all, 1, max))

write.table(res.all, file=paste("output/out", alpha.tau, alpha.mu, sigma, confounding, axis.aligned, k.tau, k.mu, p, n, "results.csv", sep="-"), row.names=FALSE, col.names=FALSE, sep=",")

}}}}}}
