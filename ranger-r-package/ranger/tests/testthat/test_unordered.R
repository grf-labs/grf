library(ranger)
library(survival)
context("ranger_unordered")

test_that("Old parameters still work", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  rf.false <- ranger(y ~ ., data = dt, num.trees = 5, write.forest = TRUE, 
                     respect.unordered.factors = FALSE)
  rf.true <- ranger(y ~ ., data = dt, num.trees = 5, write.forest = TRUE, 
                    respect.unordered.factors = TRUE)
  
  expect_null(rf.false$forest$covariate.levels)
  expect_equal(length(rf.true$forest$covariate.levels), 1)
})

test_that("If respect.unordered.factors='partition', regard characters as unordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'partition')
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'partition')
  
  expect_equal(rf.char$prediction.error, rf.fac$prediction.error)
})

test_that("If respect.unordered.factors='ignore', regard characters as ordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'ignore')
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'ignore')
  
  expect_equal(rf.char$prediction.error, rf.fac$prediction.error)
})

test_that("Error if other value for respect.unordered.factors", {
  expect_error(ranger(y ~ ., iris, num.trees = 5, respect.unordered.factors = NULL))
})

test_that("Same results if no unordered factors", {
  set.seed(100)
  rf1 <- ranger(Species ~ ., iris, num.trees = 5, respect.unordered.factors = 'ignore')
  set.seed(100)
  rf2 <- ranger(Species ~ ., iris, num.trees = 5, respect.unordered.factors = 'order')
  set.seed(100)
  rf3 <- ranger(Species ~ ., iris, num.trees = 5, respect.unordered.factors = 'partition')
  
  expect_equal(rf1$confusion.matrix, 
               rf2$confusion.matrix)
  expect_equal(rf1$confusion.matrix, 
               rf3$confusion.matrix)
})

test_that("Error if too many factors in 'partition' mode", {
  n <- 100
  dt <- data.frame(x = factor(1:100, ordered = FALSE),  
                   y = rbinom(n, 1, 0.5))
  
  expect_error(ranger(y ~ ., data = dt, num.trees = 5, respect.unordered.factors = 'partition'))
})

test_that("Survival forest with 'order' mode works", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, 
               write.forest = TRUE, respect.unordered.factors = 'order')
  expect_equal(sort(rf$forest$covariate.levels$celltype), 
               sort(levels(veteran$celltype)))
    
  pred <- predict(rf, veteran)
  expect_is(pred, "ranger.prediction")
})

test_that("maxstat splitting not working with unordered factors", {
  expect_error(ranger(Sepal.Length ~ ., iris, splitrule = "maxstat", respect.unordered.factors = "partition"))
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat", respect.unordered.factors = "partition"))
})

test_that("C splitting not working with unordered factors", {
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "C", respect.unordered.factors = "partition"))
})
