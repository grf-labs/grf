library(grf)

set.seed(1234)

test_that("Boosted regression forest improves performance vs regular forest", {

  n = 500; p = 6
  X = matrix(runif(n * p), n, p)
  mu = 2 * X[,1] * X[,2] + 3 * X[,3] + 4 * X[,4]
  Y = mu + rnorm(n)
  forest.regular <- regression_forest(X,Y)
  forest.boost <- boosted_regression_forest(X,Y)

  forest.Yhat <- predict(forest.regular)$predictions
  boost.Yhat <- predict(forest.boost)$predictions

  mse.forest <- mean((Y-forest.Yhat)^2
  mse.boost <- mean((Y- boost.Yhat))^2
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
