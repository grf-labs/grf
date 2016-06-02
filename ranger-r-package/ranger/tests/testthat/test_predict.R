## Tests for predictions

library(ranger)
context("ranger_pred")

## Tests
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

test_that("Prediction works correctly if dependent variable is not first or last", {
  dat <- iris[, c(1:2, 5, 3:4)]
  rf <- ranger(Species ~ ., dat, num.trees = 5, write.forest = TRUE)
  expect_gte(mean(predictions(predict(rf, dat)) == dat$Species), 0.9)
})


