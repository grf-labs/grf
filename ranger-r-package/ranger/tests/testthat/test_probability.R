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
  expect_is(prob$predictions, "matrix")
  expect_equal(nrow(prob$predictions), nrow(iris.test))
  expect_equal(ncol(prob$predictions), length(rg.prob$forest$levels))
})

test_that("probability estimations are between 0 and 1 and sum to 1", {
  expect_true(all(prob$predictions > -1e-5 & prob$predictions <= 1 + 1e-5))
  expect_equal(rowSums(prob$predictions), rep(1, nrow(prob$predictions)))
})

test_that("save.memory option works for probability", {
  rf <- ranger(Species ~ ., data = iris, probability = TRUE, save.memory = TRUE)
  expect_equal(rf$treetype, "Probability estimation")
})

test_that("predict works for single observations, probability prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE, probability = TRUE)
  pred <- predict(rf, head(iris, 1))
  expect_equal(names(which.max(pred$predictions)), as.character(iris[1,"Species"]))
})

test_that("Probability estimation works correctly if labels are reversed", {
  ## Simulate data
  n <- 50
  a1 <- c(rnorm(n, 3, sd = 2), rnorm(n, 8, sd = 2))
  a2 <- c(rnorm(n, 8, sd = 2), rnorm(n, 3, sd = 2))
  
  ## create labels for data
  labels <- as.factor(c(rep("0", n), rep("1", n)))
  dat <- data.frame(label = labels, a1, a2)
  
  labels.rev <- as.factor(c(rep("1", n), rep("0", n))) 
  dat.rev <- data.frame(label = labels.rev, a1, a2)
  
  ## Train
  rf <- ranger(dependent.variable.name = "label", data = dat, probability = TRUE, num.trees = 5)
  rf.rev <- ranger(dependent.variable.name = "label", data = dat.rev, probability = TRUE, num.trees = 5)
  
  ## Check predictions
  expect_gte(mean(rf$predictions[1:n, "0"]), 0.5)
  expect_gte(mean(rf$predictions[(n+1):(2*n), "1"]), 0.5)
  
  expect_gte(mean(rf.rev$predictions[1:n, "1"]), 0.5)
  expect_gte(mean(rf.rev$predictions[(n+1):(2*n), "0"]), 0.5)
})
