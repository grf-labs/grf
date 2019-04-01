library(grf)
library(testthat)
set.seed(1234)

test_that("non-vector Y throws error", {
  Y = matrix()
  expect_error(validate_Y(Y))
})

test_that("non-vector W throws error", {
  W = matrix()
  expect_error(validate_W(W))
})

test_that("factor Y throws error", {
  Y = c(1,2,3)
  Y = as.factor(Y)
  expect_error(validate_Y(Y))
})

test_that("character Y throws error", {
  Y = c(1,2,3)
  Y = as.character(Y)
  expect_error(validate_Y(Y))
})

test_that("factor W throws error", {
  W = c(0,1)
  W = as.factor(W)
  expect_error(validate_W(W))
})

test_that("character W throws error", {
  W = c(0,1)
  W = as.character(W)
  expect_error(validate_W(W))
})

test_that("W must contain only two unique values", {
  W = c(0,1,2)
  expect_error(validate_W(W))
})

test_that("X with NA throws error", {
  X = matrix(c(NA,2), nrow=1, ncol=2)
  expect_error(validate_X(X))
})

test_that("Y with NA throws error", {
  Y = c(NA,1,2)
  expect_error(validate_Y(Y))
})

test_that("W with NA throws error", {
  W = c(NA,1)
  expect_error(validate_W(W))
})

test_that("length(W) != nrow(X) throws error", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  W = c(0,1,1)
  expect_error(validate_WX(W,X))
})

test_that("length(Y) != nrow(X) throws error", {
  X = matrix(c(1,1), nrow=1, ncol=2)
  Y = c(0,1,1)
  expect_error(validate_YX(Y,X))
})
