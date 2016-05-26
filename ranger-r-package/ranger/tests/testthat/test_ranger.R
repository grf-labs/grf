library(ranger)
library(survival)
context("ranger")

test_that("no warning if data.frame has two classes", {
  dat <- iris
  class(dat) <- c("data.frame", "data.table")
  expect_that(ranger(Species ~ ., data = dat, verbose = FALSE), 
              not(gives_warning()))
})

test_that("Error if sample fraction is 0 or >1", {
  expect_that(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 0), 
              throws_error())
  expect_that(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 1.1), 
              throws_error())
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