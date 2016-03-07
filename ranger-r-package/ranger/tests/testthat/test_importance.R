library(ranger)

context("importance")

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
  expect_that(importance_pvalues(rf_p100, method = "janitza"), gives_warning())
})

test_that("Importance p-values Janitza: returns correct dimensions", {
  expect_that(
    vimp <- importance_pvalues(rf_p100, method = "janitza"), 
    gives_warning()
  )
  expect_that(vimp, is_a("matrix"))
  expect_that(dim(vimp), equals(c(104, 4)))
})

test_that("Importance p-values Janitza: error if no importance", {
  rf_none <- ranger(Species ~., iris, num.trees = 10, importance = "none", write.forest = TRUE)
  expect_that(importance_pvalues(rf_none, method = "janitza"), throws_error())
})
  
test_that("Importance p-values Janitza: error if Gini importance", {
  rf_imp <- ranger(Species ~., iris, num.trees = 10, importance = "impurity", write.forest = TRUE)
  expect_that(importance_pvalues(rf_imp, method = "janitza"), throws_error())
})

test_that("Importance p-values Janitza: error if no unimportant variables", {
  expect_that(
    expect_that(importance_pvalues(rf_p0, method = "janitza"), throws_error()), 
    gives_warning()
  )
})

test_that("Importance p-values Janitza: warning for regression", {
  rf <- ranger(Sepal.Length ~., dat_n100, num.trees = 10, importance = "permutation", write.forest = TRUE)
  expect_that(importance_pvalues(rf, method = "janitza"), gives_warning())
})

test_that("Importance p-values Janitza-Holdout: returns correct dimensions", {
  expect_that(
    vimp <- importance_pvalues(holdout_p100, method = "janitza"), gives_warning()
  )
  expect_that(vimp, is_a("matrix"))
  expect_that(dim(vimp), equals(c(104, 4)))
})

## Altmann
test_that("Importance p-values Altmann: returns correct dimensions", {
  vimp <- importance_pvalues(rf_p0, method = "altmann", formula = Species ~ ., data = iris)
  expect_that(vimp, is_a("matrix"))
  expect_that(dim(vimp), equals(c(4, 4)))
})

test_that("Importance p-values Altmann: error if no importance", {
  rf_none <- ranger(Species ~., iris, num.trees = 10, importance = "none", write.forest = TRUE)
  expect_that(importance_pvalues(rf_none, method = "altmann", formula = Species ~ ., data = iris), throws_error())
})

test_that("Importance p-values Altmann: not working for holdoutRF", {
  expect_that(importance_pvalues(holdout_p0, method = "altmann", formula = Species ~ ., data = iris), throws_error())
})

## Hold-out RF
test_that("HoldoutRF working", {
  expect_that(holdout_p0, is_a("holdoutRF"))
})

test_that("HoldoutRF ... argument working", {
  rf <- holdoutRF(Species ~., iris, num.trees = 10)
  expect_that(rf$rf1$num.trees, equals(10))
})