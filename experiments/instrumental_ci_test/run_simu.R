rm(list = ls())

setwd("~/git/grf/experiments/instrumental_ci_test2")

library(grf)

#
#
# PROBLEM CONFIGURATION
#
#

alpha.tau = 1
alpha.mu.vals = 1
sigma = 1
confounding.vals = 1
axis.aligned.vals = c(TRUE, FALSE)

k.tau.vals = c(2, 4)
k.mu = 2
p.vals = c(6, 12, 18)

n.vals = c(2000, 4000, 8000)
n.test = 1000

REPS = 20

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

forest.orthog = instrumental_forest(X, Y, W, Z, precompute.nuisance = TRUE, num.trees = 4000, ci.group.size = 10)
tau.forest.orthog = predict(forest.orthog, newdata = X.test, estimate.variance = TRUE)
tau.hat = tau.forest.orthog$predictions
sigma.hat = sqrt(tau.forest.orthog$variance.estimates)

covered = (abs(tau.hat - tau.true) / sigma.hat <= 1.96)

print(mean(covered))
return(mean(covered))

})

print(paste("CONFIG", alpha.tau, alpha.mu, sigma, confounding, axis.aligned, k.tau, k.mu, p, n, sep="-"))
print(mean(res.all))

write.table(res.all, file=paste("output/out", alpha.tau, alpha.mu, sigma, confounding, axis.aligned, k.tau, k.mu, p, n, "results.csv", sep="-"), row.names=FALSE, col.names=FALSE, sep=",")

}}}}}}
