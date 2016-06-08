## Tests for random forests for classification

library(ranger)
context("ranger_class")

## Initialize the random forest for classification
dat <- data.matrix(iris)

rg.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE)
rg.mat   <- ranger(dependent.variable.name = "Species", data = dat, write.forest = TRUE, classification = TRUE)

## Basic tests (for all random forests equal)
test_that("classification result is of class ranger with 14 elements", {
  expect_is(rg.class, "ranger")
  expect_equal(length(rg.class), 14)
})

test_that("results have 500 trees", {
  expect_equal(rg.class$num.trees, 500)
})

test_that("results have right number of independent variables", {
  expect_equal(rg.class$num.independent.variables, ncol(iris) - 1)
})

test_that("Alternative interface works for classification", {
  rf <- ranger(dependent.variable.name = "Species", data = iris)
  expect_equal(rf$treetype, "Classification")
})

test_that("Matrix interface works for classification", {
  expect_equal(rg.mat$treetype, "Classification")
  expect_equal(rg.mat$forest$independent.variable.names, colnames(iris)[1:4])
})

test_that("Matrix interface prediction works for classification", {
    expect_silent(predict(rg.mat, dat))
})

test_that("save.memory option works for classification", {
  rf <- ranger(Species ~ ., data = iris, save.memory = TRUE)
  expect_equal(rf$treetype, "Classification")
})

test_that("predict.all for classification returns numeric matrix of size trees x n", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_is(pred$predictions, "matrix")
  expect_equal(dim(pred$predictions), 
              c(nrow(iris), rf$num.trees))
})

test_that("Majority vote of predict.all for classification is equal to forest prediction", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  ## Majority vote
  pred_num <- apply(pred_trees$predictions, 1, function(x) {
    which(tabulate(x) == max(tabulate(x)))[1]
  })
  pred <- factor(pred_num, levels = 1:length(rf$forest$levels),
                 labels = rf$forest$levels)
  expect_equal(pred, pred_forest$predictions)
})

test_that("Alternative interface classification prediction works if only independent variable given, one independent variable", {
  n <- 50
  
  dt <- data.frame(x = runif(n), y = factor(rbinom(n, 1, 0.5)))
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, write.forest = TRUE)
  expect_silent(predict(rf, dt))
  expect_silent(predict(rf, dt[, 1, drop = FALSE]))
  
  dt2 <- data.frame(y = factor(rbinom(n, 1, 0.5)), x = runif(n))
  rf <- ranger(dependent.variable.name = "y", data = dt2, num.trees = 5, write.forest = TRUE)
  expect_silent(predict(rf, dt2))
  expect_silent(predict(rf, dt2[, 2, drop = FALSE]))
})

test_that("Alternative interface classification prediction works if only independent variable given, two independent variables", {
  n <- 50
  
  dt <- data.frame(x1 = runif(n), x2 = runif(n), y = factor(rbinom(n, 1, 0.5)))
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, write.forest = TRUE)
  expect_silent(predict(rf, dt))
  expect_silent(predict(rf, dt[, 1:2]))
  
  dt2 <- data.frame(y = factor(rbinom(n, 1, 0.5)), x1 = runif(n), x2 = runif(n))
  rf <- ranger(dependent.variable.name = "y", data = dt2, num.trees = 5, write.forest = TRUE)
  expect_silent(predict(rf, dt2))
  expect_silent(predict(rf, dt2[, 2:3]))
})

## Special tests for random forests for classification
test_that("predict works for single observations, classification", {
  pred <- predict(rg.class, head(iris, 1))
  expect_equal(pred$predictions, iris[1,"Species"])
})

test_that("confusion matrix is of right dimension", {
  expect_equal(dim(rg.class$confusion.matrix), 
               rep(nlevels(iris$Species), 2))
})

test_that("confusion matrix rows are the true classes", {
  expect_equal(as.numeric(rowSums(rg.class$confusion.matrix)), 
               as.numeric(table(iris$Species)))
})
