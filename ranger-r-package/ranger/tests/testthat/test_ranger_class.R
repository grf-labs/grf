##This Skript provides the tests for Rndom Forests for Classification

library(ranger)
library(survival)

context("ranger")

rg.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE)

test_that("classification result is of class ranger with 14 elements", {
  expect_that(rg.class, is_a("ranger"))
  expect_that(length(rg.class), equals(14))
})

test_that("results have 500 trees", {
  expect_that(rg.class$num.trees, equals(500))
})

test_that("results have right number of independent variables", {
  expect_that(rg.class$num.independent.variables, equals(ncol(iris) - 1))
})
