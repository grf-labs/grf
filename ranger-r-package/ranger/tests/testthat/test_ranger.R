library(ranger)
library(survival)
context("ranger")

test_that("no warning if data.frame has two classes", {
  dat <- iris
  class(dat) <- c("data.frame", "data.table")
  expect_that(ranger(Species ~ ., data = dat, verbose = FALSE), 
              not(gives_warning()))
})

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

test_that("Inbag count matrix if of right size, with replacement", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Inbag count matrix if of right size, without replacement", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, replace = FALSE, keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Inbag count matrix if of right size, with replacement, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, case.weights = runif(nrow(iris)), keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Inbag count matrix if of right size, without replacement, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, replace = FALSE, case.weights = runif(nrow(iris)), keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Error if sample fraction is 0 or >1", {
  expect_that(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 0), 
              throws_error())
  expect_that(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 1.1), 
              throws_error())
})

test_that("Number of samples is right sample fraction, replace=FALSE, default", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_that(sample.fraction, is_more_than(0.6))
  expect_that(sample.fraction, is_less_than(0.7))
})

test_that("Number of samples is right sample fraction, replace=FALSE, 0.3", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE, sample.fraction = 0.3)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_that(sample.fraction, is_more_than(0.25))
  expect_that(sample.fraction, is_less_than(0.35))
})

test_that("Number of samples is right sample fraction, replace=TRUE, default", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-1)
  
  expect_that(sample.fraction, is_more_than(expected.sample.fraction-0.05))
  expect_that(sample.fraction, is_less_than(expected.sample.fraction+0.05))
})

test_that("Number of samples is right sample fraction, replace=TRUE, 0.5", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE, sample.fraction = 0.5)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-0.5)
  
  expect_that(sample.fraction, is_more_than(expected.sample.fraction-0.05))
  expect_that(sample.fraction, is_less_than(expected.sample.fraction+0.05))
})

test_that("Number of samples is right sample fraction, replace=FALSE, 0.3, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE, sample.fraction = 0.3, case.weights = runif(nrow(iris)))
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_that(sample.fraction, is_more_than(0.25))
  expect_that(sample.fraction, is_less_than(0.35))
})

test_that("Number of samples is right sample fraction, replace=TRUE, 0.5, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE, sample.fraction = 0.5, case.weights = runif(nrow(iris)))
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-0.5)
  
  expect_that(sample.fraction, is_more_than(expected.sample.fraction-0.05))
  expect_that(sample.fraction, is_less_than(expected.sample.fraction+0.05))
})

test_that("Alternative interface classification prediction works if only independent variable given, one independent variable", {
  n <- 50
  
  dt <- data.frame(x = runif(n), y = factor(rbinom(n, 1, 0.5)))
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt), 
              not(throws_error()))
  expect_that(predict(rf, dt[, 1, drop = FALSE]), 
              not(throws_error()))
  
  dt2 <- data.frame(y = factor(rbinom(n, 1, 0.5)), x = runif(n))
  rf <- ranger(dependent.variable.name = "y", data = dt2, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt2), 
              not(throws_error()))
  expect_that(predict(rf, dt2[, 2, drop = FALSE]), 
              not(throws_error()))
})

test_that("Alternative interface classification prediction works if only independent variable given, two independent variables", {
  n <- 50
  
  dt <- data.frame(x1 = runif(n), x2 = runif(n), y = factor(rbinom(n, 1, 0.5)))
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt), 
              not(throws_error()))
  expect_that(predict(rf, dt[, 1:2]), 
              not(throws_error()))
  
  dt2 <- data.frame(y = factor(rbinom(n, 1, 0.5)), x1 = runif(n), x2 = runif(n))
  rf <- ranger(dependent.variable.name = "y", data = dt2, num.trees = 5, write.forest = TRUE)
  expect_that(predict(rf, dt2), 
              not(throws_error()))
  expect_that(predict(rf, dt2[, 2:3]), 
              not(throws_error()))
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

test_that("No error if survival tree without OOB observations", {
  dat <- data.frame(time = c(1,2), status = c(0,1), x = c(1,2))
  expect_that(ranger(Surv(time, status) ~ ., dat, num.trees = 1, num.threads = 1), 
              not(throws_error()))
})

test_that("as.factor() in formula works", {
  n <- 20
  dt <- data.frame(x = runif(n), y = rbinom(n, 1, 0.5))
  expect_that(ranger(as.factor(y) ~ ., data = dt, num.trees = 5, write.forest = TRUE), 
              not(throws_error()))
})

test_that("If respect.unordered.factors=TRUE, regard characters as unordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = TRUE)
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = TRUE)
  
  expect_that(rf.char$prediction.error, equals(rf.fac$prediction.error))
})

test_that("If respect.unordered.factors=FALSE, regard characters as ordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = FALSE)
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = FALSE)
  
  expect_that(rf.char$prediction.error, equals(rf.fac$prediction.error))
})