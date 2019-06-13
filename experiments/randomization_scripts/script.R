library(grf)

rfseed <- 1234
set.seed(rfseed)

# Generate data for different forests
n = 200; p = 10
num.trees = 200
num.threads = 2
X = matrix(rnorm(n*p), n, p)
X.test = matrix(0, 101, p)
X.test[,1] = seq(-2, 2, length.out = 101)
Z = (runif(n) > 0.5) * rlnorm(n, sdlog = 2)
W = rbinom(n, 1, 0.4 + 0.2 * (X[,1] > 0))
Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n) + Z

# Regression forest
rf <- regression_forest(X, Y, num.trees = num.trees, seed = rfseed, num.threads=num.threads)
print(head(predict(rf)$predictions, 5))

# Causal forest
cf = causal_forest(X, Y, W, seed = rfseed, num.threads=num.threads)
print(head(predict(cf)$predictions, 5))

# ATE
ow_ate <- average_treatment_effect(cf, target.sample = "overlap")
print(ow_ate)
