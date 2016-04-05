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

test_that("predict.all for classification returns numeric matrix of size trees x n", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_that(pred$predictions, is_a("matrix"))
  expect_that(dim(pred$predictions), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Majority vote of predict.all for classification is equal to forest prediction", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  ## Majority vote
  pred_num <- apply(pred_trees$predictions, 1, function(x) {
    which(tabulate(x) == max(tabulate(x)))
  })
  pred <- factor(pred_num, levels = 1:length(rf$forest$levels),
                 labels = rf$forest$levels)
  expect_that(pred, equals(pred_forest$predictions))
})

##Special tests for random forests for classification
test_that("predict works for single observations, classification", {
  pred <- predict(rg.class, head(iris, 1))
  expect_that(pred$predictions, equals(iris[1,"Species"]))
})

test_that("confusion matrix is of right dimension", {
  expect_that(dim(rg.class$confusion.matrix), 
              equals(rep(nlevels(iris$Species), 2)))
})

test_that("confusion matrix rows are the true classes", {
  expect_that(as.numeric(rowSums(rg.class$confusion.matrix)), 
              equals(as.numeric(table(iris$Species))))
})