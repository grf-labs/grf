library(grf)

set.seed(1234)

test_that("local linear prediction gives reasonable estimates", {
  f <- function(x) {
    x[1] + 2 * x[2] + 2 * x[3]**2
  }
  n <- 600
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  MU <- apply(X, FUN = f, MARGIN = 1)
  Y <- MU + rnorm(n)

  forest <- regression_forest(X, Y, num.trees = 500)
  preds.grf.oob <- predict(forest)
  preds.ll.oob <- predict(forest, linear.correction.variables = 1:p, ll.lambda = 0)

  mse.grf.oob <- mean((preds.grf.oob$predictions - MU)^2)
  mse.ll.oob <- mean((preds.ll.oob$predictions - MU)^2)

  expect_lt(mse.ll.oob, 1)
  expect_lt(mse.ll.oob, mse.grf.oob / 2)

  X.test <- matrix(rnorm(n * p), n, p)
  MU.test <- apply(X.test, FUN = f, MARGIN = 1)

  preds.grf <- predict(forest, X.test)
  preds.ll <- predict(forest, X.test, linear.correction.variables = 1:p, ll.lambda = 0.1)

  mse.grf <- mean((preds.grf$predictions - MU.test)^2)
  mse.ll <- mean((preds.ll$predictions - MU.test)^2)

  expect_lt(mse.ll, 1)
  expect_lt(mse.ll, mse.grf / 1.5)
})

test_that("linear correction variables function as expected", {
  f <- function(x) {
    x[1] + 2 * x[2] + 2 * x[3]**2
  }
  n <- 400
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  MU <- apply(X, FUN = f, MARGIN = 1)
  Y <- MU + rnorm(n)

  forest <- regression_forest(X, Y, num.trees = 500)
  preds <- predict(forest, linear.correction.variables = 1:20)
  mse <- mean((preds$predictions - MU)^2)

  preds.selected <- predict(forest, linear.correction.variables = 1:3)
  mse.selected <- mean((preds.selected$predictions - MU)^2)

  expect_lt(mse.selected, mse / 1.5)
})

test_that("local linear forest tuning returns lambda and decreases prediction error", {
  n <- 400
  p <- 5
  sigma <- 5

  mu <- function(x) {
    log(1 + exp(6 * x[1]))
  }

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  truth <- apply(X, FUN = mu, MARGIN = 1)
  Y <- truth + sigma * rnorm(n)

  forest <- ll_regression_forest(X, Y)

  lambda <- tune_ll_regression_forest(forest)$lambda.min
  expect_true(is.numeric(lambda))
  expect_true(length(lambda) == 1)

  preds.tuned <- predict(forest, linear.correction.variables = 1:p, ll.lambda = lambda)$predictions
  mse.tuned <- mean((preds.tuned - truth)^2)

  preds.untuned <- predict(forest, linear.correction.variables = 1:p, ll.lambda = 0.1)$predictions
  mse.untuned <- mean((preds.untuned - truth)^2)

  expect_lt(mse.tuned, 0.75 * mse.untuned)
})

test_that("default local linear forest predict and regression forest predict with local.linear = TRUE are the same", {
  n <- 400
  p <- 5
  sigma <- 1

  mu <- function(x) {
    log(1 + exp(6 * x[1]))
  }

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  truth <- apply(X, FUN = mu, MARGIN = 1)
  Y <- truth + sigma * rnorm(n)

  forest <- regression_forest(X, Y)
  preds <- predict(forest, linear.correction.variables = 1:5, lambda = 0.1)$predictions

  ll.forest <- ll_regression_forest(X, Y)
  ll.preds <- predict(ll.forest, ll.lambda = 0.1)$predictions

  average.difference <- mean((ll.preds - preds)**2)

  expect_lt(average.difference, 0.05)
})

test_that("local linear predict returns local linear predictions even without tuning parameters", {
  n <- 50
  p <- 5
  sigma <- 1

  mu <- function(x) {
    log(1 + exp(6 * x[1]))
  }

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  truth <- apply(X, FUN = mu, MARGIN = 1)
  Y <- truth + sigma * rnorm(n)

  forest <- ll_regression_forest(X, Y, num.trees = 50)
  preds <- predict(forest)

  ll.indicator <- !is.null(preds$ll.lambda)
  expect_true(ll.indicator)

  forest <- ll_regression_forest(X, Y, num.trees = 50, enable.ll.split = TRUE)
  preds <- predict(forest)

  ll.indicator <- !is.null(preds$ll.lambda)
  expect_true(ll.indicator)
})

test_that("local linear confidence intervals have reasonable coverage", {
  mu <- function(x) {
    log(1 + exp(6 * x))
  }

  n <- 400
  p <- 10
  sigma <- sqrt(20)

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  truth <- mu(X[, 1])
  Y <- truth + sigma * rnorm(n)

  forest <- ll_regression_forest(X, Y, num.trees = 500, ci.group.size = 5)
  preds <- predict(forest, linear.correction.variables = 1, ll.lambda = 1, estimate.variance = TRUE)

  expect_true(all(preds$variance.estimates > 0))

  df <- data.frame(
    predictions = preds$predictions,
    upper = preds$predictions + 1.96 * sqrt(preds$variance.estimates),
    lower = preds$predictions - 1.96 * sqrt(preds$variance.estimates)
  )

  percent_llf <- mean(df$lower <= truth & truth <= df$upper)
  expect_gt(percent_llf, 0.8)
})


test_that("local linear confidence intervals match regression forest with large lambda", {
  mu <- function(x) {
    log(1 + exp(6 * x))
  }

  n <- 80
  p <- 4
  sigma <- sqrt(20)

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  truth <- mu(X[, 1])
  Y <- truth + sigma * rnorm(n)

  forest <- regression_forest(X, Y, num.trees = 80, ci.group.size = 2)

  preds.rf <- predict(forest, estimate.variance = TRUE)
  preds.llf <- predict(forest,
    linear.correction.variables = 1,
    ll.lambda = 10000000, estimate.variance = TRUE
  )

  expect_lt(max(abs(preds.llf$predictions - preds.rf$predictions)), 10^-5)
  expect_lt(max(abs(preds.llf$variance.estimates - preds.rf$variance.estimates)), 10^-5)
})

test_that("local linear predictions are correct without noise", {
  n <- 80
  p <- 2

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  mu <- rowSums(X)
  Y <- mu

  forest <- regression_forest(X, Y, num.trees = 80, ci.group.size = 2)

  preds.rf <- predict(forest)$predictions
  preds.llf <- predict(forest, linear.correction.variables = 1:p, ll.lambda = 0)$predictions

  expect_lt(mean((preds.llf - mu)^2), 10^-10)
  expect_gt(mean((preds.rf - mu)^2), 10^-2)

  forest <- ll_regression_forest(X, Y, num.trees = 80, enable.ll.split = TRUE, ci.group.size = 2)
  preds.llf.splits <- predict(forest, linear.correction.variables = 1:p, ll.lambda = 0)$predictions
  expect_lt(mean((preds.llf.splits - mu)^2), 10^-10)
})

test_that("prediction with and without CIs are the same", {
  n <- 200
  p <- 4

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  mu <- 0.9 * exp(X[, 1])
  Y <- mu + rnorm(n)

  forest <- ll_regression_forest(X, Y, num.trees = 800, ci.group.size = 2)

  preds.rf <- predict(forest, ll.lambda = 1, estimate.variance = FALSE)$predictions
  preds.rf2 <- predict(forest, ll.lambda = 1, estimate.variance = TRUE)$predictions
  expect_lt(max(abs(preds.rf - preds.rf2)), 10^-10)
})

test_that("output of tune local linear forest is consistent with prediction output", {
  n <- 200
  p <- 4

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  mu <- 0.9 * exp(X[, 1])
  Y <- mu + rnorm(n)

  forest <- ll_regression_forest(X, Y, num.trees = 400)

  tune.out <- tune_ll_regression_forest(forest)

  ll.min <- tune.out$lambdas[1]
  pred.ll.min <- predict(forest, ll.lambda = ll.min)$predictions
  expect_lt(max(abs(tune.out$oob.predictions[, 1] - pred.ll.min)), 10^-6)

  ll.max <- tune.out$lambdas[length(tune.out$lambdas)]
  pred.ll.max <- predict(forest, ll.lambda = ll.max)$predictions
  expect_lt(max(abs(tune.out$oob.predictions[, length(tune.out$lambdas)] - pred.ll.max)), 10^-6)
})

test_that("local linear forests with local linear splits are numeric and include variance estimates", {
  n <- 200
  p <- 4

  X <- matrix(runif(n * p, -1, 1), nrow = n)
  X.test <- matrix(runif(n * p, -1, 1), nrow = n)
  mu <- 0.9 * exp(X[, 1])
  Y <- mu + rnorm(n)

  forest <- ll_regression_forest(X, Y, num.trees = 800, enable.ll.split = TRUE, ci.group.size = 2)
  preds.oob <- predict(forest, estimate.variance = TRUE)
  preds.test <- predict(forest, X.test, estimate.variance = TRUE)

  expect_true(is.numeric(preds.oob$predictions))
  expect_true(is.numeric(preds.test$predictions))

  expect_equal(n, length(preds.oob$variance.estimates))
  expect_equal(n, length(preds.test$variance.estimates))
})

test_that("local linear splits reduce early splits on linear trends", {
   n <- 200
   p <- 6

   X <- matrix(runif(n * p, -1, 1), nrow = n)
   mu <- 0.9 * exp(X[, 1]) + 2 * X[,6]
   Y <- mu + rnorm(n)

   ll.forest <- ll_regression_forest(X, Y, enable.ll.split = TRUE)
   ll.split.freq <- split_frequencies(ll.forest, 1)

   forest = regression_forest(X, Y)
   split.freq <- split_frequencies(forest, 1)

   expect_lt(split.freq[1,1], ll.split.freq[1,1] / 3)
   expect_gt(split.freq[1,6], ll.split.freq[1,6] * 3)
})

test_that("local linear splits improve predictions in a simple case", {
   n <- 600
   p <- 5

   X <- matrix(rnorm(n * p, 0, 1), nrow = n)
   MU <- X[,1] + X[,2] + X[,3]*X[,4]
   Y <- MU + rnorm(n)

   forest <- regression_forest(X, Y, num.trees = 500)
   preds.grf.splits.oob <- predict(forest, linear.correction.variables = 1:p, ll.lambda = 0.1)

   ll.forest <- ll_regression_forest(X, Y, num.trees = 500, enable.ll.split = TRUE)
   preds.ll.splits.oob <- predict(ll.forest, linear.correction.variables = 1:p, ll.lambda = 0.1)

   mse.grf.splits.oob <- mean((preds.grf.splits.oob$predictions - MU)^2)
   mse.ll.splits.oob <- mean((preds.ll.splits.oob$predictions - MU)^2)

   expect_lt(mse.ll.splits.oob / mse.grf.splits.oob, 0.9)
})

test_that("local linear split regulating works in a simple case", {
   n <- 2000
   p <- 5

   X <- matrix(runif(n * p, -1, 1), nrow = n)
   mu <- 0.9 * exp(X[, 1]) + 2 * X[,2] + 2 * X[,3] * X[,4]
   Y <- mu + rnorm(n)

   forest.regulate <- ll_regression_forest(X, Y, num.trees = 500, enable.ll.split = TRUE)
   forest <- ll_regression_forest(X, Y, num.trees = 500, enable.ll.split = TRUE, ll.split.cutoff = 0)

   preds.regulate <- predict(forest.regulate)$predictions
   preds <- predict(forest)$predictions

   mse.preds.regulate <- mean((preds.regulate - mu)^2)
   mse.preds <- mean((preds - mu)^2)
   expect_lt(mse.preds.regulate / mse.preds, 0.5)
})
