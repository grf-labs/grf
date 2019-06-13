library(grf)

seed <- 1000
set.seed(seed)

test_that("Boosted regression forest improves performance vs regular forest", {
  seed <- 1000
  set.seed(seed)
  n <- 200;
  p <- 6
  X <- matrix(runif(n * p), n, p)
  mu <- 2 * X[,1] * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y <- mu + rnorm(n)
  forest.regular <- regression_forest(X, Y , seed = seed)
  forest.boost <- boosted_regression_forest(X, Y, seed = seed)

  forest.Yhat <- predict(forest.regular)$predictions
  boost.Yhat <- predict(forest.boost)$predictions

  mse.forest <- mean((Y - forest.Yhat)^2)
  mse.boost <- mean((Y - boost.Yhat)^2)
  expect_true(mse.boost < mse.forest)
})

test_that("Boosted forest takes user specified number of steps",  {
  seed <- 1000
  set.seed(seed)
  n <- 100; p <- 6
  X <- matrix(runif(n * p), n, p)
  mu <- 2 * X[,1] * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y <- mu + rnorm(n)
  forest.boost <- boosted_regression_forest(X,Y,boost.steps=2, seed = seed)
  expect_equal(2,length(forest.boost$forests))
})

test_that("Under cross-validation, errors decrease each step", {
  n<-200; p<-6
  X <- matrix(runif(n * p), n, p)
  mu <- 2 * X[,1]^2 * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y <- mu + rnorm(n)
  forest.boost <- boosted_regression_forest(X, Y, seed = seed)
  errors <- unlist(forest.boost$error)
  if(length(errors)>1) {
    improves <- errors[2:length(errors)] - errors[1:(length(errors)-1)]
    expect_true(all(improves<0))
  }
})

test_that("boost.error.reduction validation works", {
  seed <- 1000
  set.seed(seed)
  n<-200; p<-6
  X <- matrix(runif(n * p), n, p)
  mu <- 2 * X[,1]^2 * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y <- mu + rnorm(n)
  expect_error(forest.boost <- boosted_regression_forest(X,Y,boost.error.reduction=1.5))
})

test_that("OOB prediction is close to actual out of sample error", {
  seed <- 1000
  set.seed(seed)
  n <- 5000;
  p <- 6
  X <- matrix(runif(n * p), n, p)
  mu <- 2 * X[,1]^2 * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y <- mu + rnorm(n)
  train <- 1:4000
  test <- 4000:5000
  forest.boost <- boosted_regression_forest(X[train,], Y[train], seed = seed)
  OOB.error <- mean((forest.boost$predictions - Y[train])^2)

  test.pred <- predict(forest.boost, X[test,])$predictions
  test.error <- mean((test.pred - Y[test])^2)

  expect_true(abs(test.error - OOB.error) < 0.15)
})
