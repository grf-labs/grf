library(ranger)

context("importance")

## GenABEL data
dat_gwaa <- readRDS("../test_gwaa.Rds")

## 0 noise variables
rf_p0 <- ranger(Species ~., iris, num.trees = 100, 
                importance = "permutation", write.forest = TRUE)
holdout_p0 <- holdoutRF(Species ~., iris, num.trees = 10)

## 100 noise variables
n <- nrow(iris)
p <- 100
noise <- replicate(p, rnorm(n))
colnames(noise) <- paste0("noise", 1:p)
dat_n100 <- cbind(iris, noise)

rf_p100 <- ranger(Species ~., dat_n100, num.trees = 100,
                  importance = "permutation", write.forest = TRUE)
holdout_p100 <- holdoutRF(Species ~., dat_n100, num.trees = 100)

## Janitza
test_that("Importance p-values Janitza: warning if few negative importance values", {
  expect_warning(importance_pvalues(rf_p100, method = "janitza"))
})

test_that("Importance p-values Janitza: returns correct dimensions", {
  expect_warning(vimp <- importance_pvalues(rf_p100, method = "janitza"))
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(104, 4))
})

test_that("Importance p-values Janitza: error if no importance", {
  rf_none <- ranger(Species ~., iris, num.trees = 10, importance = "none", write.forest = TRUE)
  expect_error(importance_pvalues(rf_none, method = "janitza"))
})
  
test_that("Importance p-values Janitza: error if Gini importance", {
  rf_imp <- ranger(Species ~., iris, num.trees = 10, importance = "impurity", write.forest = TRUE)
  expect_error(importance_pvalues(rf_imp, method = "janitza"))
})

test_that("Importance p-values Janitza: error if no unimportant variables", {
  expect_warning(expect_error(
    importance_pvalues(rf_p0, method = "janitza")))
})

test_that("Importance p-values Janitza: warning for regression", {
  rf <- ranger(Sepal.Length ~., dat_n100, num.trees = 10, importance = "permutation", write.forest = TRUE)
  expect_warning(importance_pvalues(rf, method = "janitza"))
})

test_that("Importance p-values Janitza-Holdout: returns correct dimensions", {
  expect_warning(vimp <- importance_pvalues(holdout_p100, method = "janitza"))
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(104, 4))
})

## Altmann
test_that("Importance p-values Altmann: returns correct dimensions", {
  vimp <- importance_pvalues(rf_p0, method = "altmann", formula = Species ~ ., data = iris)
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(4, 4))
})

test_that("Importance p-values Altmann: error if no importance", {
  rf_none <- ranger(Species ~., iris, num.trees = 10, importance = "none", write.forest = TRUE)
  expect_error(importance_pvalues(rf_none, method = "altmann", formula = Species ~ ., data = iris))
})

test_that("Importance p-values Altmann: not working for holdoutRF", {
  expect_error(importance_pvalues(holdout_p0, method = "altmann", formula = Species ~ ., data = iris))
})

## Hold-out RF
test_that("HoldoutRF working", {
  expect_is(holdout_p0, "holdoutRF")
})

test_that("HoldoutRF working with GenABEL data", {
  holdout_gwaa <- holdoutRF(CHD ~., dat_gwaa, num.trees = 10)
  expect_is(holdout_p0, "holdoutRF")
})

test_that("HoldoutRF ... argument working", {
  rf <- holdoutRF(Species ~., iris, num.trees = 10)
  expect_equal(rf$rf1$num.trees, 10)
})
