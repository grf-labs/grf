library(grf)

library(testthat)
set.seed(1234)

test_that("non-vector observation throws error", {
  X <- matrix(c(1, 1), nrow = 1, ncol = 2)
  Y <- matrix()
  expect_error(validate_observations(Y, X))
})

test_that("factor Y throws error", {
  X <- matrix(c(1, 1), nrow = 1, ncol = 2)
  Y <- c(1, 2, 3)
  Y <- as.factor(Y)
  expect_error(validate_observations(Y, X))
})

test_that("character Y throws error", {
  X <- matrix(c(1, 1), nrow = 1, ncol = 2)
  Y <- c(1, 2, 3)
  Y <- as.character(Y)
  expect_error(validate_observations(Y, X))
})

test_that("X with NA throws error", {
  X <- matrix(c(NA, 2), nrow = 1, ncol = 2)
  expect_error(validate_X(X))
})

test_that("Y with NA throws error", {
  X <- matrix(c(1, 1), nrow = 1, ncol = 2)
  Y <- c(NA, 1, 2)
  expect_error(validate_observations(Y, X))
})

test_that("validate_observations processes list of vectors", {
  X <- matrix(c(1, 1), nrow = 1, ncol = 2)
  Y <- c(1)
  W <- c(2, 3)
  expect_error(validate_observations(list(Y, W), X))
})

test_that("length(Y) != nrow(X) throws error", {
  X <- matrix(c(1, 1), nrow = 1, ncol = 2)
  Y <- c(0, 1, 1)
  expect_error(validate_observations(Y, X))
})

test_that("an observation matrix with 1 column is accepted", {
  X <- matrix(c(0, 0, 1, 1, 2, 2), nrow = 3, ncol = 2)
  Y <- matrix(c(0, 1, 2), nrow = 3, ncol = 1)
  colnames(Y) <- "outcome"
  Y <- validate_observations(Y, X)
  expect_true(is.vector(Y))
})

test_that("create_train_matrices handles data.frame, matrix, and NULL inputs equally", {
  Xm <- matrix(rnorm(100), 20, 5)
  Xd <- as.data.frame(Xm)
  Y <- matrix(rnorm(20))

  data1_d <- create_train_matrices(Xd)
  data1_m <- create_train_matrices(Xm)
  data2_d <- create_train_matrices(Xd, Y)
  data2_m <- create_train_matrices(Xm, Y)
  data3_d <- create_train_matrices(Xd, Y, NULL)
  data3_m <- create_train_matrices(Xm, Y, NULL)

  # Checking for equality of elements
  # (note expect_equal does not work here)
  expect_true(all(data1_d$default == data1_m$default))
  expect_true(all(data2_d$default == data2_m$default))
  expect_true(all(data3_d$default == data3_m$default))
})

test_that("providing sample.weights when equalize.cluster.weights is TRUE is not accepted", {
  equalize.cluster.weights <- TRUE
  clusters = c(1, 1, 2, 2)
  sample.weights = c(2, 2, 1, 1)

  expect_error(validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights))
})
