##This skript provides the tests for character vector in data

library(ranger)
library(survival)
context("ranger")

##Initialize random forests
dat <- iris
dat$Test <- paste0("AA",as.character(1:nrow(dat)))

##Tests
test_that("no warning if character vector in data", {
  expect_that(ranger(Species ~ ., data = dat, verbose = FALSE), 
              not(gives_warning()))
})

test_that("no error if character vector in data, prediction", {
  rf <- ranger(Species~., dat, write.forest = TRUE)
  expect_that(predict(rf, dat),
              not(throws_error()))
})

test_that("no warning if character vector in data, alternative interface", {
  expect_that(ranger(dependent.variable.name = "Species", data = dat, verbose = FALSE), 
              not(gives_warning()))
})

test_that("no error if character vector in data, alternative interface, prediction", {
  rf <- ranger(dependent.variable.name = "Species", data = dat, verbose = FALSE, write.forest = TRUE)
  expect_that(predict(rf, dat),
              not(throws_error()))
})
