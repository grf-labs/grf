library(ranger)
library(survival)

context("ranger")

rg.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE)
rg.reg <- ranger(Sepal.Length ~ ., data = iris, verbose = FALSE, write.forest = TRUE)
rg.surv <- ranger(Surv(time, status) ~ ., data = veteran, verbose = FALSE, write.forest = TRUE)

test_that("classification result is of class ranger with 14 elements", {
  expect_that(rg.class, is_a("ranger"))
  expect_that(length(rg.class), equals(14))
})

test_that("regression result is of class ranger with 14 elements", {
  expect_that(rg.reg, is_a("ranger"))
  expect_that(length(rg.reg), equals(14))
})

test_that("survival result is of class ranger with 15 elements", {
  expect_that(rg.surv, is_a("ranger"))
  expect_that(length(rg.surv), equals(15))
})

test_that("results have 500 trees", {
  expect_that(rg.class$num.trees, equals(500))
  expect_that(rg.reg$num.trees, equals(500))
  expect_that(rg.surv$num.trees, equals(500))
})

test_that("results have right number of independent variables", {
  expect_that(rg.class$num.independent.variables, equals(ncol(iris) - 1))
  expect_that(rg.reg$num.independent.variables, equals(ncol(iris) - 1))
  expect_that(rg.surv$num.independent.variables, equals(ncol(veteran) - 2))
})

test_that("unique death times in survival result is right", {
  expect_that(rg.surv$unique.death.times, equals(sort(unique(veteran$time))))
})

test_that("importance measures work", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "impurity")
  expect_that(rg.imp$variable.importance, is_a("numeric"))
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation")
  expect_that(rg.imp$variable.importance, is_a("numeric"))
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation", scale.permutation.importance = TRUE)
  expect_that(rg.imp$variable.importance, is_a("numeric"))
})

test_that("gini importance is larger than 1", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "impurity")
  expect_that(rg.imp$variable.importance[1], is_more_than(1))
})

test_that("unscaled importance is smaller than 1", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation", scale.permutation.importance = FALSE)
  expect_that(rg.imp$variable.importance[1], is_less_than(1))
})

test_that("scaled importance is larger than 1", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation", scale.permutation.importance = TRUE)
  expect_that(rg.imp$variable.importance[1], is_more_than(1))
})

test_that("probability estimations are a matrix with correct size", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]

  rg.prob <- ranger(Species ~ ., data = iris.train, write.forest = TRUE, probability = TRUE)
  prob <- predict(rg.prob, iris.test)

  expect_that(prob$predictions, is_a("matrix"))
  expect_that(nrow(prob$predictions), equals(nrow(iris.test)))
  expect_that(ncol(prob$predictions), equals(length(rg.prob$forest$levels)))
})

test_that("probability estimations are between 0 and 1 and sum to 1", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]

  rg.prob <- ranger(Species ~ ., data = iris.train, write.forest = TRUE, probability = TRUE)
  prob <- predict(rg.prob, iris.test)

  expect_that(all(prob$predictions > -1e-5 & prob$predictions <= 1 + 1e-5), is_true())
  expect_that(rowSums(prob$predictions), equals(rep(1, nrow(prob$predictions))))
})

test_that("predict returns good prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE)
  pred <- predict(rf, iris)
  expect_that(mean(iris$Species == predictions(pred)), is_more_than(0.9))
})
