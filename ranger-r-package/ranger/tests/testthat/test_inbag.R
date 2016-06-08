## Tests for inbag functions

library(ranger)
context("ranger_inbag")

## Tests
test_that("Inbag count matrix if of right size, with replacement", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE)
  expect_equal(dim(data.frame(rf$inbag.counts)), 
              c(nrow(iris), rf$num.trees))
})

test_that("Inbag count matrix if of right size, without replacement", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, replace = FALSE, keep.inbag = TRUE)
  expect_equal(dim(data.frame(rf$inbag.counts)), 
              c(nrow(iris), rf$num.trees))
})

test_that("Inbag count matrix if of right size, with replacement, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, case.weights = runif(nrow(iris)), keep.inbag = TRUE)
  expect_equal(dim(data.frame(rf$inbag.counts)), 
              c(nrow(iris), rf$num.trees))
})

test_that("Inbag count matrix if of right size, without replacement, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, replace = FALSE, case.weights = runif(nrow(iris)), keep.inbag = TRUE)
  expect_equal(dim(data.frame(rf$inbag.counts)), 
              c(nrow(iris), rf$num.trees))
})


test_that("Number of samples is right sample fraction, replace=FALSE, default", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_gt(sample.fraction, 0.6)
  expect_lt(sample.fraction, 0.7)
})

test_that("Number of samples is right sample fraction, replace=FALSE, 0.3", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE, sample.fraction = 0.3)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_gt(sample.fraction, 0.25)
  expect_lt(sample.fraction, 0.35)
})

test_that("Number of samples is right sample fraction, replace=TRUE, default", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-1)
  
  expect_gt(sample.fraction, expected.sample.fraction-0.05)
  expect_lt(sample.fraction, expected.sample.fraction+0.05)
})

test_that("Number of samples is right sample fraction, replace=TRUE, 0.5", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE, sample.fraction = 0.5)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-0.5)
  
  expect_gt(sample.fraction, expected.sample.fraction-0.05)
  expect_lt(sample.fraction, expected.sample.fraction+0.05)
})

test_that("Number of samples is right sample fraction, replace=FALSE, 0.3, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE, sample.fraction = 0.3, case.weights = runif(nrow(iris)))
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_gt(sample.fraction, 0.25)
  expect_lt(sample.fraction, 0.35)
})

test_that("Number of samples is right sample fraction, replace=TRUE, 0.5, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE, sample.fraction = 0.5, case.weights = runif(nrow(iris)))
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-0.5)
  
  expect_gt(sample.fraction, expected.sample.fraction-0.05)
  expect_lt(sample.fraction, expected.sample.fraction+0.05)
})
