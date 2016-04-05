##This skript provides the tests for random forests for regression

library(ranger)
library(survival)
context("ranger")

##Initialize the random forest for regression
rg.reg <- ranger(Sepal.Length ~ ., data = iris, verbose = FALSE, write.forest = TRUE)

##Basic tests (for all random forests equal)
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
##Special tests for random forests for regression