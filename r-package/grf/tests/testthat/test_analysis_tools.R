library(grf)

set.seed(1234)

test_that("examining a tree gives reasonable results", {
  p <- 40
  n <- 500

  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))

  qrf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry = p, min.node.size = 10, sample.fraction = 0.632)
  quantile.tree <- get_tree(qrf, 500)

  num.nodes <- length(quantile.tree$nodes)
  expect_lt(num.nodes, n)

  split.vars <- unlist(sapply(quantile.tree$nodes, function(node) node$split_variable))
  expect_true(all(split.vars >= 0 && split.vars <= 40))

  left.children <- unlist(sapply(quantile.tree$nodes, function(node) node$left_child))
  expect_equal(seq(2, num.nodes, 2), left.children)
  right.children <- unlist(sapply(quantile.tree$nodes, function(node) node$right_child))
  expect_equal(seq(3, num.nodes, 2), right.children)
})

test_that("leaf samples are indexed correctly", {
  p <- 40
  n <- 500

  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))

  qrf <- quantile_forest(X, Y, sample.fraction = 1.0, honesty = FALSE)
  quantile.tree <- get_tree(qrf, 1)

  samples <- unlist(sapply(quantile.tree$nodes, function(node) node$samples))
  expect_equal(seq(1, n), sort(samples))
})

test_that("draw samples are indexed correctly", {
  p <- 40
  n <- 500

  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))

  forest <- regression_forest(X, Y)
  forest.tree <- get_tree(forest, 1)

  leaf_nodes <- Filter(f = function(x) x$is_leaf, forest.tree$nodes)

  # This should contain all in-bag data
  estimation_and_split_sample <- forest.tree$drawn_samples

  # This is the estimation sample. It should be contained in the vector above
  estimation_sample <- unlist(Map(f = function(x) x$samples, leaf_nodes))

  # This shouldn't contain anything...
  should_be_empty <- setdiff(estimation_sample, estimation_and_split_sample)
  expect_equal(length(should_be_empty), 0)
})




test_that("tree indexes are one-indexed", {
  p <- 40
  n <- 500

  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))

  qrf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry = p, min.node.size = 10, sample.fraction = 0.632)

  expect_error(get_tree(qrf, 0))
  get_tree(qrf, n)
})

test_that("calculating variable importance gives reasonable results", {
  p <- 40
  n <- 500

  i <- 5
  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- rnorm(n) * (1 + (X[, i] > 0))

  qrf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), mtry = p, min.node.size = 10, sample.fraction = 0.632)

  var.importance <- variable_importance(qrf)

  expect_equal(length(var.importance), p)
  expect_equal(which.max(var.importance), i)
})

test_that("computing sample weights gives reasonable results", {
  p <- 10
  n <- 100

  X <- matrix(2 * runif(n * p) - 1, n, p)
  Y <- (X[, 1] > 0) + 2 * rnorm(n)

  rrf <- regression_forest(X, Y, mtry = p)

  sample.weights.oob <- get_sample_weights(rrf)
  expect_equal(nrow(sample.weights.oob), n)
  expect_equal(ncol(sample.weights.oob), n)

  row.sums.oob <- apply(sample.weights.oob, 1, sum)
  expect_true(all(row.sums.oob - 1.0 < 1e-10))

  n.test <- 103
  X.test <- matrix(2 * runif(n.test * p) - 1, n.test, p)

  sample.weights <- get_sample_weights(rrf, X.test)
  expect_equal(nrow(sample.weights), n.test)
  expect_equal(ncol(sample.weights), n)

  row.sums <- apply(sample.weights, 1, sum)
  expect_true(all(row.sums - 1.0 < 1e-10))
})

test_that("regression forest leaf nodes contains 'avg Y' only", {
  n <- 50
  p <- 1
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * rnorm(n)
  r.forest <- regression_forest(X, Y)
  r.tree <- get_tree(r.forest,1)
  for(n in r.tree$nodes){
    if(n$is_leaf){
      expect_false(is.null(names(n$leaf_stats)))
      expect_setequal(names(n$leaf_stats), c("avg Y"))
    }
  }
})

test_that("causal forest leaf nodes contains 'avg Y' and 'avg W' only", {
  p <- 4
  n <- 100
  X <- matrix(runif(n * p), n, p)
  Y <- runif(n)
  W <- rbinom(n, 1, 0.5)
  c.forest <- causal_forest(X, Y, W, num.trees = 100)
  c.tree <- get_tree(c.forest,1)
  for(n in c.tree$nodes){
    if(n$is_leaf){
      expect_false(is.null(names(n$leaf_stats)))
      expect_setequal(names(n$leaf_stats), c("avg Y", "avg W"))
    }
  }
})

test_that("instrumental forest leaf nodes contains 'avg Y', 'avg W', and 'avg Z' only", {
  p <- 6
  n <- 200

  X <- matrix(rnorm(n * p), n, p)
  
  eps <- rnorm(n)
  Z <- rbinom(n, 1, 2 / 3)
  filter <- rbinom(n, 1, 1 / (1 + exp(-1 * eps)))
  W <- Z * filter
  
  tau <- apply(X[, 1:2], 1, function(xx) sum(pmax(0, xx)))
  mu <- apply(X[, 2 + 1:2], 1, function(xx) sum(pmax(0, xx)))
  
  Y <- (2 * W - 1) / 2 * tau + mu + eps
  
  iv.forest <- instrumental_forest(X, Y, W, Z, num.trees = 100)
  iv.tree <- get_tree(iv.forest,1)

  for(n in iv.tree$nodes){
    if(n$is_leaf){
      expect_false(is.null(names(n$leaf_stats)))
      expect_setequal(names(n$leaf_stats), c("avg Y", "avg W", "avg Z"))
    }
  }
})

