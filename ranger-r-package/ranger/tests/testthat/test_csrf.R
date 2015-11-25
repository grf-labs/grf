library(ranger)
library(survival)

context("case-specific RF")

test_that("csrf classification returns predictions", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]
  
  pred <- csrf(Species ~ ., training_data = iris.train, test_data = iris.test, 
               params1 = list(num.trees = 50), 
               params2 = list(num.trees = 10))
  
  expect_that(pred, is_a("factor"))
  expect_that(length(pred), equals(50))
})

test_that("csrf regression returns predictions", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]
  
  pred <- csrf(Sepal.Length ~ ., training_data = iris.train, test_data = iris.test, 
               params1 = list(num.trees = 50), 
               params2 = list(num.trees = 10))
  
  expect_that(pred, is_a("numeric"))
  expect_that(length(pred), equals(50))
})