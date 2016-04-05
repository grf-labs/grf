##This skript provides the tests for predictions

library(ranger)
library(survival)
context("ranger")


test_that("predict returns good prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE)
  pred <- predict(rf, iris)
  expect_that(mean(iris$Species == predictions(pred)), is_more_than(0.9))
})

test_that("case weights work", {
  expect_that(ranger(Species ~ ., iris, num.trees = 5, case.weights = rep(1, nrow(iris))), 
              not(throws_error()))
  ## Should only predict setosa now
  weights <- c(rep(1, 50), rep(0, 100))
  rf <- ranger(Species ~ ., iris, num.trees = 5, case.weights = weights, write.forest = TRUE)
  pred <- predict(rf, iris)$predictions
  expect_that(all(pred == "setosa"), is_true())
})


