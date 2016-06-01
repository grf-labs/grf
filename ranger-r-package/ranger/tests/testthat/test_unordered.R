library(ranger)
library(survival)
context("ranger_unordered")

test_that("If unordered.factors='partition', regard characters as unordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, unordered.factors='partition')
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, unordered.factors='partition')
  
  expect_that(rf.char$prediction.error, equals(rf.fac$prediction.error))
})

test_that("If unordered.factors='ignore', regard characters as ordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, unordered.factors='ignore')
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, unordered.factors='ignore')
  
  expect_that(rf.char$prediction.error, equals(rf.fac$prediction.error))
})
