library(ranger)

context("importance")

## 0 noise variables
rf_p0 <- ranger(Species ~., iris, num.trees = 100, 
                importance = "permutation", write.forest = TRUE)

## 100 noise variables
n <- nrow(iris)
p <- 100
noise <- replicate(p, rnorm(n))
colnames(noise) <- paste0("noise", 1:p)
dat_n100 <- cbind(iris, noise)

rf_p100 <- ranger(Species ~., dat_n100, num.trees = 100,
                  importance = "permutation", write.forest = TRUE)

## 1000 noise variables
n <- nrow(iris)
p <- 1000
noise <- replicate(p, rnorm(n))
colnames(noise) <- paste0("noise", 1:p)
dat_n1000 <- cbind(iris, noise)

rf_p1000 <- ranger(Species ~., dat_n1000, num.trees = 1000,
                   importance = "permutation", write.forest = TRUE)

## Janitza
test_that("Importance p-values Janitza: warning if few negative importance values", {
  expect_that(importance_pvalues(rf_p100, method = "janitza"), gives_warning())
})

test_that("Importance p-values Janitza: returns correct dimensions", {
  vimp <- importance_pvalues(rf_p1000, method = "janitza")
  expect_that(vimp, is_a("matrix"))
  expect_that(dim(vimp), equals(c(1004, 2)))
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
  expect_that(importance_pvalues(rf_p0, method = "janitza"), throws_error())
})

test_that("Importance p-values Janitza: warning for regression", {
  rf <- ranger(Sepal.Length ~., dat_n100, num.trees = 10, importance = "permutation", write.forest = TRUE)
  expect_that(importance_pvalues(rf, method = "janitza"), gives_warning())
})

## Altmann
test_that("Importance p-values Altmann: returns correct dimensions", {
  vimp <- importance_pvalues(rf_p0, method = "altmann", formula = Species ~ ., data = iris)
  expect_that(vimp, is_a("matrix"))
  expect_that(dim(vimp), equals(c(4, 2)))
})

test_that("Importance p-values Altmann: error if no importance", {
  rf_none <- ranger(Species ~., iris, num.trees = 10, importance = "none", write.forest = TRUE)
  expect_that(importance_pvalues(rf_none, method = "altmann", formula = Species ~ ., data = iris), throws_error())
})