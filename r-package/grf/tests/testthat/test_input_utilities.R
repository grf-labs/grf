library(grf)
library(testthat)
set.seed(1234)

test_that("non-vector observation throws error", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  Y = matrix()
  expect_error(validate_observations(Y,X))
})

test_that("factor Y throws error", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  Y = c(1,2,3)
  Y = as.factor(Y)
  expect_error(validate_observations(Y,X))
})

test_that("character Y throws error", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  Y = c(1,2,3)
  Y = as.character(Y)
  expect_error(validate_observations(Y,X))
})

test_that("X with NA throws error", {
  X = matrix(c(NA,2), nrow=1, ncol=2)
  expect_error(validate_X(X))
})

test_that("Y with NA throws error", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  Y = c(NA,1,2)
  expect_error(validate_observations(Y,X))
})

test_that("validate_observations processes list of vectors", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  Y = c(1)
  W = c(2,3)
  expect_error(validate_observations(list(Y,W),X))
})

test_that("length(Y) != nrow(X) throws error", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  Y = c(0,1,1)
  expect_error(validate_observations(Y,X))
})

test_that("an observation matrix with 1 column is accepted", {
  X = matrix(c(0, 0, 1, 1, 2, 2), nrow=3, ncol=2)
  Y = matrix(c(0,1,2), nrow=3, ncol=1)
  colnames(Y) = "outcome"

  Y = validate_observations(Y, X)
  expect_true(is.vector(Y))
})
