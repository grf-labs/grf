## Tests for split select weights

library(ranger)
context("ranger_splitweights")

## Tests
test_that("split select weights work", {
  expect_that(ranger(Species ~ ., iris, num.trees = 5, split.select.weights = c(0.1, 0.2, 0.3, 0.4)), 
              not(throws_error()))
  expect_that(ranger(Species ~ ., iris, num.trees = 5, split.select.weights = c(0.1, 0.2, 0.3)), 
              throws_error())
})

test_that("Tree-wise split select weights work", {
  num.trees <- 5
  weights <- replicate(num.trees, runif(ncol(iris)-1), simplify = FALSE)
  expect_that(ranger(Species ~ ., iris, num.trees = num.trees, split.select.weights = weights), 
              not(throws_error()))
  
  weights <- replicate(num.trees+1, runif(ncol(iris)-1), simplify = FALSE)
  expect_that(ranger(Species ~ ., iris, num.trees = num.trees, split.select.weights = weights), 
              throws_error())
})
