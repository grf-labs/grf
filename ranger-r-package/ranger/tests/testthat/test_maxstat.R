library(ranger)
library(survival)
context("ranger_maxstat")

test_that("maxstat splitting works for survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.4)
})

test_that("maxstat splitting works for regression", {
  rf <- ranger(Sepal.Length ~ ., iris, splitrule = "maxstat")
  expect_is(rf, "ranger")
  expect_gt(rf$r.squared, 0.5)
})

test_that("maxstat splitting, alpha out of range throws error", {
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat", alpha = -1))
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat", alpha = 2))
})