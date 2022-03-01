# grf weighting benchmark script
# this script benchmarks statistical performance
# for sample weighted causal forest and rank average treatment effect
library(grf)
set.seed(1)

n <- 40000
p <- 4
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)
U <- rbinom(n, 1, 0.75)
TAU <- U * pmax(X[, 1], 0)
Y <- TAU * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

G <- 1 - 2/3 * U
cc <- which(rbinom(n, 1, G) == 1)

tau.forest.1 <- causal_forest(X, Y, W, W.hat = 0.5)
tau.forest.2 <- causal_forest(X[cc,], Y[cc], W[cc], W.hat = 0.5)
tau.forest.3 <- causal_forest(X[cc,], Y[cc], W[cc], sample.weights = 1/G[cc], W.hat = 0.5)

prio <- X[,1]
rate1 <- rank_average_treatment_effect(tau.forest.1, prio)
rate2 <- rank_average_treatment_effect(tau.forest.2, prio[cc])
rate3 <- rank_average_treatment_effect(tau.forest.3, prio[cc])

rate1
rate2
rate3

# These should be different
c(rate1$estimate, rate2$estimate)

# These should be similar
c(rate1$estimate, rate3$estimate)

# These should be similar (because in cc sample have 1/2 instead of 3/4 samples with non-zero TAU)
c(rate1$estimate * 2 / 3, rate2$estimate)
