##This skript provides the tests for random forests for classification

library(ranger)
library(survival)
context("ranger")

##Initialize the random forest for classification
rg.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE)

##Basic tests (for all random forests equal)
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

test_that("Alternative interface works for classification", {
  rf <- ranger(dependent.variable.name = "Species", data = iris)
  expect_that(rf$treetype, equals("Classification"))
})

test_that("Matrix interface works for classification", {
  rf <- ranger(dependent.variable.name = "Species", data = data.matrix(iris), write.forest = TRUE, classification = TRUE)
  expect_that(rf$treetype, equals("Classification"))
  expect_that(rf$forest$independent.variable.names, equals(colnames(iris)[1:4]))
})

test_that("Matrix interface prediction works for classification", {
  dat <- data.matrix(iris)
  rf <- ranger(dependent.variable.name = "Species", data = dat, write.forest = TRUE, classification = TRUE)
  expect_that(predict(rf, dat), not(throws_error()))
})

test_that("save.memory option works for classification", {
  rf <- ranger(Species ~ ., data = iris, save.memory = TRUE)
  expect_that(rf$treetype, equals("Classification"))
})
##Special tests for random forests for classification
test_that("predict works for single observations, classification", {
  pred <- predict(rg.class, head(iris, 1))
  expect_that(pred$predictions, equals(iris[1,"Species"]))
})