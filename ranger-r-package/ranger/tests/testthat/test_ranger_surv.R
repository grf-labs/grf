## Tests for random forests for survival analysis

library(ranger)
library(survival)
context("ranger_surv")

## Initialize the random forest for survival analysis
rg.surv <- ranger(Surv(time, status) ~ ., data = veteran, verbose = FALSE, write.forest = TRUE)

## Basic tests (for all random forests equal)
test_that("survival result is of class ranger with 16 elements", {
  expect_that(rg.surv, is_a("ranger"))
  expect_that(length(rg.surv), equals(16))
})

test_that("results have 500 trees", {
  expect_that(rg.surv$num.trees, equals(500))
})

test_that("results have right number of independent variables", {
  expect_that(rg.surv$num.independent.variables, equals(ncol(veteran) - 2))
})

test_that("Alternative interface works for survival", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran)
  expect_that(rf$treetype, equals("Survival"))
})

test_that("Matrix interface works for survival", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = data.matrix(veteran), write.forest = TRUE)
  expect_that(rf$treetype, equals("Survival"))
  expect_that(rf$forest$independent.variable.names, equals(colnames(veteran)[c(1:2, 5:8)]))
})

test_that("Matrix interface prediction works for survival", {
  dat <- data.matrix(veteran)
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = dat, write.forest = TRUE)
  expect_that(predict(rf, dat), not(throws_error()))
})

test_that("predict works for single observations, survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, write.forest = TRUE)
  pred <- predict(rf, head(veteran, 1))
  expect_that(length(pred$survival), equals(length(rf$unique.death.times)))
})

## Special tests for random forests for survival analysis
test_that("unique death times in survival result is right", {
  expect_that(rg.surv$unique.death.times, equals(sort(unique(veteran$time))))
})

test_that("C-index splitting works", {
  rf <- ranger(Surv(time, status) ~ ., data = veteran, verbose = FALSE, 
               splitrule = "C")
  expect_that(rf$treetype, equals("Survival"))
})

test_that("No error if survival tree without OOB observations", {
  dat <- data.frame(time = c(1,2), status = c(0,1), x = c(1,2))
  expect_that(ranger(Surv(time, status) ~ ., dat, num.trees = 1, num.threads = 1), 
              not(throws_error()))
})