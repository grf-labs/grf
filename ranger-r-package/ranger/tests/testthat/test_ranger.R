library(ranger)
library(survival)
context("ranger")

## GenABEL
if (!requireNamespace("GenABEL", quietly = TRUE)) {
  stop("Package GenABEL is required for testing ranger completely. Please install it.", call. = FALSE)
} else {
  dat.gwaa <- readRDS("../test_gwaa.Rds")
  rg.gwaa <- ranger(CHD ~ ., data = dat.gwaa, verbose = FALSE, write.forest = TRUE)
}

test_that("classification gwaa rf is of class ranger with 14 elements", {
  expect_is(rg.gwaa, "ranger")
  expect_equal(length(rg.gwaa), 14)
})

test_that("Matrix interface works for Probability estimation", {
  rf <- ranger(dependent.variable.name = "Species", data = data.matrix(iris), write.forest = TRUE, probability = TRUE)
  expect_equal(rf$treetype, "Probability estimation")
  expect_equal(rf$forest$independent.variable.names, colnames(iris)[1:4])
})

test_that("Matrix interface prediction works for Probability estimation", {
  dat <- data.matrix(iris)
  rf <- ranger(dependent.variable.name = "Species", data = dat, write.forest = TRUE, probability = TRUE)
  expect_silent(predict(rf, dat))
})

test_that("no warning if data.frame has two classes", {
  dat <- iris
  class(dat) <- c("data.frame", "data.table")
  expect_silent(ranger(Species ~ ., data = dat, verbose = FALSE))
})

test_that("Error if sample fraction is 0 or >1", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 0))
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 1.1))
})

test_that("as.factor() in formula works", {
  n <- 20
  dt <- data.frame(x = runif(n), y = rbinom(n, 1, 0.5))
  expect_silent(ranger(as.factor(y) ~ ., data = dt, num.trees = 5, write.forest = TRUE))
})

test_that("maxstat splitting works for survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat")
  expect_is(rf, "ranger")
})

test_that("maxstat splitting, alpha out of range throws error", {
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat", alpha = -1))
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat", alpha = 2))
})

test_that("holdout mode holding out data with 0 weight", {
  weights <- rbinom(nrow(iris), 1, 0.5)
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, sample.fraction = 0.632*mean(weights), 
               holdout = TRUE, keep.inbag = TRUE)
  inbag <- data.frame(rf$inbag.counts)
  expect_true(all(inbag[weights == 0, ] == 0))
})

test_that("holdout mode uses holdout OOB data", {
  weights <- rbinom(nrow(iris), 1, 0.5)
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, sample.fraction = 0.632*mean(weights), 
               holdout = TRUE, keep.inbag = TRUE)
  expect_false(any(is.na(rf$predictions[weights == 0])))
  expect_true(all(is.na(rf$predictions[weights == 1])))
})

test_that("holdout mode not working if no weights", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, importance = "permutation", holdout = TRUE))
})

test_that("holdout mode: no OOB prediction if no 0 weights", {
  weights <- runif(nrow(iris))
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, 
               holdout = TRUE, keep.inbag = TRUE)
  expect_true(all(is.na(rf$predictions)))
})

test_that("Probability estimation works for empty classes", {
  expect_silent(rf <- ranger(Species ~., iris[1:100,],  num.trees = 5, probability = TRUE))
})

