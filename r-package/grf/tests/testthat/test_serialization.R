library(grf)

save_and_load = function(forest) {
  save('forest', file='forest.Rdata')
  load('forest.Rdata')
  forest
}

identical_predictions = function(forest.a, forest.b, X.new=NULL) {
  if(is.null(X.new)) { expect_true(all(predict(forest.a) == predict(forest.b))) }
  else { expect_true(all(predict(forest.a, X.new) == predict(forest.b, X.new))) }
}

test_that("regression forests make the same predictions before and after save/load.", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[,1] + 0.1 * rnorm(n)
  X.new <- matrix(rnorm(n * p), n, p)

  forest = regression_forest(X, Y, serialize=TRUE)
  identical_predictions(forest, save_and_load(forest))
  identical_predictions(forest, save_and_load(forest), X.new)
})

test_that("local linear forests make the same predictions before and after save/load.", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[,1] + 0.1 * rnorm(n)
  X.new <- matrix(rnorm(n * p), n, p)

  forest = local_linear_forest(X, Y, serialize=TRUE)
  identical_predictions(forest, save_and_load(forest))
  identical_predictions(forest, save_and_load(forest), X.new)
})

test_that("quantile forests make the same predictions before and after save/load", {
  p = 10
  n = 500
  i = 5
  X = matrix(2 * runif(n * p) - 1, n, p)
  Y = rnorm(n)  + 100 * (X[,i] > 0)
  X.new = matrix(2 * runif(n * p) - 1, n, p)

  forest = quantile_forest(X, Y)
  identical_predictions(forest, save_and_load(forest), X.new)
})

test_that("causal forests make the same predictions before and after save/load", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- X[,1] * (2 * W - 1) + 0.1 * rnorm(n)

  forest = causal_forest(X, Y, W)
  identical_predictions(forest, save_and_load(forest))
})

test_that("instrumental forests make the same predictions before and after save/load", {
  p = 6
  n = 2000
  n.test = 2000
  X = matrix(rnorm(n * p), n, p)
  eps = rnorm(n)
  Z = rbinom(n, 1, 2/3)
  filter = rbinom(n, 1, 1/(1 + exp(-1 * eps)))
  W = Z * filter
  tau = apply(X[,1:2], 1, function(xx) sum(pmax(0, xx)))
  mu = apply(X[,2 + 1:2], 1, function(xx) sum(pmax(0, xx)))
  Y = (2 * W - 1) / 2 * tau + mu + eps

  forest = instrumental_forest(X, Y, W, Z, num.trees=4000)
  identical_predictions(forest, save_and_load(forest))
})


