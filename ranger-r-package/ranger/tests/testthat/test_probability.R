## Tests for random forests for probability estimation

library(ranger)
context("ranger_prob")

## Initialize random forest
train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
iris.train <- iris[train.idx, ]
iris.test <- iris[-train.idx, ]

rg.prob <- ranger(Species ~ ., data = iris.train, write.forest = TRUE, probability = TRUE)
prob <- predict(rg.prob, iris.test)

## Tests
test_that("probability estimations are a matrix with correct size", {
  expect_that(prob$predictions, is_a("matrix"))
  expect_that(nrow(prob$predictions), equals(nrow(iris.test)))
  expect_that(ncol(prob$predictions), equals(length(rg.prob$forest$levels)))
})

test_that("probability estimations are between 0 and 1 and sum to 1", {
  expect_that(all(prob$predictions > -1e-5 & prob$predictions <= 1 + 1e-5), is_true())
  expect_that(rowSums(prob$predictions), equals(rep(1, nrow(prob$predictions))))
})

test_that("save.memory option works for probability", {
  rf <- ranger(Species ~ ., data = iris, probability = TRUE, save.memory = TRUE)
  expect_that(rf$treetype, equals("Probability estimation"))
})

test_that("predict works for single observations, probability prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE, probability = TRUE)
  pred <- predict(rf, head(iris, 1))
  expect_that(names(which.max(pred$predictions)), equals(as.character(iris[1,"Species"])))
})
