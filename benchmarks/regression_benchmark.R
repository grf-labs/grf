library(grf)
library(microbenchmark)
library(ranger)

p = 6
n = 1000

X = matrix(2 * runif(n * p) - 1, n, p)
Y = (X[,1] > 0) + 2 * rnorm(n)
data <- data.frame(Y, X)

ticks = 101
X.test = matrix(0, ticks, p)
xvals = seq(-1, 1, length.out = ticks)
X.test[,1] = xvals
truth = xvals > 0

test.data = data.frame(X.test)

bench.train = microbenchmark(
  grf.forest = regression_forest(X, Y, num.trees = 1000),
  ranger.forest = ranger(Y ~ ., data = data, num.trees = 1000),
  times = 3, unit = "s")

grf.forest = regression_forest(X, Y, num.trees = 1000)
ranger.forest = ranger(Y ~ ., data = data, num.trees = 1000)

bench.pred = microbenchmark(
  grf.pred = predict(grf.forest, X.test),
  ranger.pred = predict(ranger.forest, test.data),
  times = 3, unit = "s")
