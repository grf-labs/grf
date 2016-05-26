## Tests for random forests for regression

library(ranger)
context("ranger_reg")

## Initialize the random forest for regression
rg.reg <- ranger(Sepal.Length ~ ., data = iris, verbose = FALSE, write.forest = TRUE)

## Basic tests (for all random forests equal)
test_that("regression result is of class ranger with 14 elements", {
  expect_that(rg.reg, is_a("ranger"))
  expect_that(length(rg.reg), equals(14))
})

test_that("results have 500 trees", {
  expect_that(rg.reg$num.trees, equals(500))
})

test_that("results have right number of independent variables", {
  expect_that(rg.reg$num.independent.variables, equals(ncol(iris) - 1))
})

test_that("Alternative interface works for regression", {
  rf <- ranger(dependent.variable.name = "Sepal.Length", data = iris)
  expect_that(rf$treetype, equals("Regression"))
})

test_that("Matrix interface works for regression", {
  rf <- ranger(dependent.variable.name = "Sepal.Length", data = data.matrix(iris), write.forest = TRUE)
  expect_that(rf$treetype, equals("Regression"))
  expect_that(rf$forest$independent.variable.names, equals(colnames(iris)[2:5]))
})

test_that("Matrix interface prediction works for regression", {
  dat <- data.matrix(iris)
  rf <- ranger(dependent.variable.name = "Sepal.Length", data = dat, write.forest = TRUE)
  expect_that(predict(rf, dat), not(throws_error()))
})

test_that("save.memory option works for regression", {
  rf <- ranger(Sepal.Length ~ ., data = iris, save.memory = TRUE)
  expect_that(rf$treetype, equals("Regression"))
})

test_that("predict.all for regression returns numeric matrix of size n x trees", {
  rf <- ranger(Petal.Width ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_that(pred$predictions, is_a("matrix"))
  expect_that(dim(pred$predictions), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Mean of predict.all for regression is equal to forest prediction", {
  rf <- ranger(Petal.Width ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  expect_that(rowMeans(pred_trees$predictions), equals(pred_forest$predictions))
})

test_that("Alternative interface regression prediction works if only independent variable given, one independent variable", {
  n <- 50
  
  dt <- data.frame(x = runif(n), y = rbinom(n, 1, 0.5))
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt), 
              not(throws_error()))
  expect_that(predict(rf, dt[, 1, drop = FALSE]), 
              not(throws_error()))
  
  dt2 <- data.frame(y = rbinom(n, 1, 0.5), x = runif(n))
  rf <- ranger(dependent.variable.name = "y", data = dt2, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt2), 
              not(throws_error()))
  expect_that(predict(rf, dt2[, 2, drop = FALSE]), 
              not(throws_error()))
})

test_that("Alternative interface regression prediction works if only independent variable given, two independent variables", {
  n <- 50
  
  dt <- data.frame(x1 = runif(n), x2 = runif(n), y = rbinom(n, 1, 0.5))
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt), 
              not(throws_error()))
  expect_that(predict(rf, dt[, 1:2]), 
              not(throws_error()))
  
  dt2 <- data.frame(y = rbinom(n, 1, 0.5), x1 = runif(n), x2 = runif(n))
  rf <- ranger(dependent.variable.name = "y", data = dt2, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt2), 
              not(throws_error()))
  expect_that(predict(rf, dt2[, 2:3]), 
              not(throws_error()))
})

test_that("Alternative interface regression prediction: Results not all the same", {
  n <- 50
  
  dt <- data.frame(x = runif(n), y = rbinom(n, 1, 0.5))
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, write.forest = TRUE)
  expect_that(diff(range(predict(rf, dt)$predictions)), is_more_than(0))
  expect_that(diff(range(predict(rf, dt[, 1, drop = FALSE])$predictions)), is_more_than(0))
  
  dt2 <- data.frame(y = rbinom(n, 1, 0.5), x = runif(n))
  rf <- ranger(dependent.variable.name = "y", data = dt2, num.trees = 5, write.forest = TRUE)
  expect_that(diff(range(predict(rf, dt2)$predictions)), is_more_than(0))
  expect_that(diff(range(predict(rf, dt2[, 2, drop = FALSE])$predictions)), is_more_than(0))
})

## Special tests for random forests for regression
