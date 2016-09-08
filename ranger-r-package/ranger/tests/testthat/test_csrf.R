library(ranger)
library(survival)

context("case-specific RF")

test_that("csrf classification returns predictions", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]
  
  pred <- csrf(Species ~ ., training_data = iris.train, test_data = iris.test, 
               params1 = list(num.trees = 10), 
               params2 = list(num.trees = 3))
  
  expect_is(pred, "factor")
  expect_equal(length(pred), nrow(iris)/3)
})

test_that("csrf regression returns predictions", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]
  
  pred <- csrf(Sepal.Length ~ ., training_data = iris.train, test_data = iris.test, 
               params1 = list(num.trees = 10), 
               params2 = list(num.trees = 3))
  
  expect_is(pred, "numeric")
  expect_equal(length(pred), nrow(iris)/3)
})
