library(grf)

set.seed(1)

test_that("Boosted regression forest improves performance vs regular forest", {

  n = 500; p = 6
  X = matrix(runif(n * p), n, p)
  mu = 2 * X[,1] * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y = mu + rnorm(n)
  forest.regular <- regression_forest(X,Y)
  forest.boost <- boosted_regression_forest(X,Y)

  forest.Yhat <- predict(forest.regular)$predictions
  boost.Yhat <- predict(forest.boost)$predictions

  mse.forest <- mean((Y-forest.Yhat)^2)
  mse.boost <- mean((Y- boost.Yhat)^2)
  expect_true(mse.boost < mse.forest)
})


test_that("Boosted forest takes user specified number of steps",  {
  n = 100; p = 6
  X = matrix(runif(n * p), n, p)
  mu = 2 * X[,1] * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y = mu + rnorm(n)
  forest.boost <- boosted_regression_forest(X,Y,num.steps=2)
  expect_equal(2,length(forest.boost$forests))
})

 test_that("Tolerance above 1 increases length of boosting", {
   n=200; p=6
   X = matrix(runif(n * p), n, p)
   mu = 2 * X[,1]^2 * X[,2] + 3 * X[,3] + 4 * X[,4]
   Y = mu + rnorm(n)
   forest.boost.bad <- boosted_regression_forest(X,Y,tolerance =1.2)
   forest.boost <- boosted_regression_forest(X,Y,num.steps=2)
   expect_true(length(forest.boost$forests) < length(forest.boost.bad$forests))
 })

 test_that("Under cross-validation, errors decrease each step", {
   n=200; p=6
   X = matrix(runif(n * p), n, p)
   mu = 2 * X[,1]^2 * X[,2] + 3 * X[,3] + 4 * X[,4]
   Y = mu + rnorm(n)
   forest.boost <- boosted_regression_forest(X,Y,num.steps=2)
   errors <- unlist(lapply(forest.boost$debiased.errors,mean))
   improves <- errors[2:length(errors)] - errors[1:(length(errors)-1)]
   expect_true(sum(improves>0)==0)

 })
