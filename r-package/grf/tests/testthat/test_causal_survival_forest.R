library(grf)

test_that("causal survival forest is well-calibrated", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, n.mc = 5000, dgp = "simple1")
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 500)
  cs.pred <- predict(cs.forest)

  mse.oob <- mean((cs.pred$predictions - data$cate)^2)

  X.test <- matrix(0.5, 10, p)
  X.test[, 1] <- seq(0, 1, length.out = 10)
  true.effect.test <- generate_causal_survival_data(10, p, Y.max = 1, X = X.test, n.mc = 5000, dgp = "simple1")$cate
  cs.pred.test <- predict(cs.forest, X.test)
  mse.test <- mean((cs.pred.test$predictions - true.effect.test)^2)

  expect_lt(mse.oob, 0.01)
  expect_lt(mse.test, 0.01)
})

test_that("causal survival forest with complete non-censored data is identical to causal forest", {
  n <- 500
  p <- 5
  X <- matrix(runif(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y.max <- 1
  failure.time <- pmin(rexp(n) * X[, 1] + W, Y.max)
  censor.time <- 999 * runif(n)
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)
  
  cs.forest <- causal_survival_forest(X, Y, W, D)
  pp.cs <- predict(cs.forest, estimate.variance = TRUE)
  cf <- causal_forest(X, Y, W)
  pp.cf <- predict(cf, estimate.variance = TRUE)
  
  expect_lt(mean((pp.cs$predictions - pp.cf$predictions)^2), 0.0005)
  expect_equal(mean(pp.cs$predictions), mean(pp.cf$predictions), tolerance = 0.01)
  expect_lt(mean((pp.cs$variance.estimates - pp.cf$variance.estimates)^2), 0.0005)
  expect_equal(mean(pp.cs$variance.estimates), mean(pp.cf$variance.estimates), tolerance = 0.01)
})

test_that("causal survival forest predictions are kernel weighted correctly", {
  # Test that the sufficient statistic approach to solving the forest weighted
  # estimating equation is internally consistent.
  n <- 250
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, n.mc = 1, dgp = "simple1")
  sample.weights <- sample(c(1, 10), n, TRUE)
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 250)
  cs.forest.weighted <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 250, sample.weights = sample.weights)
  x1 <- data$X[1, , drop = FALSE]
  cs.pred <- predict(cs.forest, x1)$predictions
  cs.pred.weighted <- predict(cs.forest.weighted, x1)$predictions
  forest.weights <- get_forest_weights(cs.forest, x1)[1, ]
  forest.weights.weighted <- get_forest_weights(cs.forest.weighted, x1)[1, ]

  theta1 <- sum(forest.weights * cs.forest$eta$numerator) /
    sum(forest.weights * cs.forest$eta$denominator)
  theta1.weighted <- sum(forest.weights.weighted * sample.weights * cs.forest.weighted$eta$numerator) /
    sum(forest.weights.weighted * sample.weights * cs.forest.weighted$eta$denominator)

  expect_equal(cs.pred, theta1)
  expect_equal(cs.pred.weighted, theta1.weighted)
})

test_that("causal survival forest variance estimates are decent", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, n.mc = 5000, dgp = "simple1")
  true.effect <- data$cate
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 500)
  cs.pred <- predict(cs.forest, estimate.variance = TRUE)
  z.score.oob <- abs(cs.pred$predictions - true.effect) / sqrt(cs.pred$variance.estimates)
  cate.coverage.oob <- mean(z.score.oob <= 1.96)
  expect_gte(cate.coverage.oob, 0.7)

  X.test <- matrix(0.5, 10, p)
  X.test[, 1] <- seq(0, 1, length.out = 10)
  true.effect.test <- generate_causal_survival_data(10, p, Y.max = 1, X = X.test, n.mc = 5000, dgp = "simple1")$cate
  cs.pred.test <- predict(cs.forest, X.test, estimate.variance = TRUE)
  z.score.test <- abs(cs.pred.test$predictions - true.effect.test) / sqrt(cs.pred.test$variance.estimates)
  cate.coverage.test <- mean(z.score.test <= 1.96)
  expect_gte(cate.coverage.test, 0.6)

  # Duplicate some samples
  sample.weights <- sample(c(1, 2), n, TRUE)
  cs.forest.weighted <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 500, sample.weights = sample.weights)
  cs.pred.weighted <- predict(cs.forest.weighted, estimate.variance = TRUE)
  z.score.weighted <- abs(cs.pred.weighted$predictions - true.effect) / sqrt(cs.pred.weighted$variance.estimates)
  cate.coverage.oob.weighted <- mean(z.score.weighted <= 1.96)
  expect_gte(cate.coverage.oob.weighted, 0.7)
})

test_that("sample weighted causal survival forest is invariant to scaling", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  failure.time <- exp(0.5 * X[, 1] + X[, 2] * W) * rexp(n)
  censor.time <- 2
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)
  weights <- runif(n)

  # Predictions are invariant to sample weight scaling
  sf <- causal_survival_forest(X, Y, W, D, num.trees = 100, sample.weights = weights, seed = 1)
  sf2 <- causal_survival_forest(X, Y, W, D, num.trees = 100, sample.weights = weights * 2, seed = 1)
  p1 <- predict(sf)
  p2 <- predict(sf2)

  expect_equal(p1$predictions, p2$predictions)

  # Variance estimates are invariant to sample weight scaling
  v1 <- predict(sf, estimate.variance = TRUE)
  v2 <- predict(sf2, estimate.variance = TRUE)

  expect_equal(v1$variance.estimates, v2$variance.estimates)
})

test_that("a causal survival forest trained on dense data is near identical to forest trained on coarser grid", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, n.mc = 1, dgp = "simple1")
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 500)
  cs.pred <- predict(cs.forest)

  Y.grid <- seq(min(data$Y), max(data$Y), length.out = 25)
  cs.forest.grid <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 500,
                                           failure.times = Y.grid)
  cs.pred.grid <- predict(cs.forest.grid)

  expect_lt(mean((cs.pred$predictions - cs.pred.grid$predictions)^2), 0.001)
})

test_that("nuisance argument handling in a causal survival forest works as expected", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, n.mc = 1, dgp = "simple1")
  cs.forest.default <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 200)
  cs.pred <- predict(cs.forest.default)

  sf.survival <- survival_forest(cbind(data$X, data$W), data$Y, data$D, num.trees = 200)
  # S1.hat and S0.hat are not proper OOB estimates but should be OK for this testing purpose.
  S1.hat <- predict(sf.survival, cbind(data$X, rep(1, n)))
  E1.hat <- expected_survival(S1.hat$predictions, sf.survival$failure.times)
  S0.hat <- predict(sf.survival, cbind(data$X, rep(0, n)))
  E0.hat <- expected_survival(S0.hat$predictions, sf.survival$failure.times)
  S.hat <- predict(sf.survival, failure.times = sort(unique(data$Y)))
  sf.censor <- survival_forest(cbind(data$X, data$W), data$Y, 1 - data$D, num.trees = 200)
  C.hat <- predict(sf.censor, failure.times = sort(unique(data$Y)), prediction.type = "Nelson-Aalen")

  cs.forest1 <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 200,
                                       E1.hat = E1.hat)
  cs.forest2 <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 200,
                                       E0.hat = E0.hat)
  cs.forest3 <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 200,
                                       S.hat = S.hat$predictions)
  cs.forest4 <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 200,
                                       C.hat = C.hat$predictions)
  cs.forest5 <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 200,
                                       E1.hat = E1.hat,
                                       E0.hat = E0.hat,
                                       S.hat = S.hat$predictions,
                                       C.hat = C.hat$predictions)

  expect_lt(mean((predict(cs.forest1)$predictions - cs.pred$predictions)^2), 0.005)
  expect_lt(mean((predict(cs.forest2)$predictions - cs.pred$predictions)^2), 0.005)
  expect_lt(mean((predict(cs.forest3)$predictions - cs.pred$predictions)^2), 0.005)
  expect_lt(mean((predict(cs.forest4)$predictions - cs.pred$predictions)^2), 0.005)
  expect_lt(mean((predict(cs.forest5)$predictions - cs.pred$predictions)^2), 0.005)
})

test_that("causal survival forest works as expected with missing values", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, n.mc = 1, dgp = "simple1")
  nmissing <- 200
  data$X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  csf <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 500)

  # MIA with data duplication
  Xl <- data$X
  Xr <- data$X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)
  csf.mia <- causal_survival_forest(X.mia, data$Y, data$W, data$D, num.trees = 500)

  mse.oob.diff <- mean((predict(csf)$predictions - predict(csf.mia)$predictions)^2)
  mse.diff <- mean((predict(csf, data$X)$predictions - predict(csf.mia, X.mia)$predictions)^2)

  expect_equal(mse.oob.diff, 0, tolerance = 0.001)
  expect_equal(mse.diff, 0, tolerance = 0.001)
})

# This characterization test locks in behavior for the default causal survival forest.
# It is done here in addition to ForestCharacterizationTest.cpp as the computation of
# nuisance components involves a fair amount of work in R.
test_that("causal survival forest has not changed ", {
  set.seed(42)
  n <- 500
  p <- 5
  dgp <- "simple1"
  data <- generate_causal_survival_data(n = n, p = p, dgp = dgp)
  cs.forest <- causal_survival_forest(round(data$X, 2), round(data$Y, 2), data$W, data$D,
                                      num.trees = 50, seed = 42, num.threads = 4)

  expected.predictions.oob <- as.numeric(readLines("data/causal_survival_oob_predictions.csv"))
  expected.predictions <- as.numeric(readLines("data/causal_survival_predictions.csv"))

  # Update with:
  # write.table(predict(cs.forest)$predictions, file = "data/causal_survival_oob_predictions.csv", row.names = FALSE, col.names = FALSE)
  # write.table(predict(cs.forest, round(data$X, 2))$predictions, file = "data/causal_survival_predictions.csv", row.names = FALSE, col.names = FALSE)


  expect_equal(predict(cs.forest)$predictions, expected.predictions.oob)
  expect_equal(predict(cs.forest, round(data$X, 2))$predictions, expected.predictions)

  # Same forest with constrained event grid.
  failure.times <- seq(min(data$Y), max(data$Y), length.out = 5)
  cs.forest.grid <- causal_survival_forest(round(data$X, 2), round(data$Y, 2), data$W, data$D,
                                           failure.times = failure.times,
                                           num.trees = 50, seed = 42, num.threads = 4)

  expected.predictions.oob.grid <- as.numeric(readLines("data/causal_survival_oob_predictions_grid.csv"))
  expected.predictions.grid <- as.numeric(readLines("data/causal_survival_predictions_grid.csv"))

  # Update with:
  # write.table(predict(cs.forest.grid)$predictions, file = "data/causal_survival_oob_predictions_grid.csv", row.names = FALSE, col.names = FALSE)
  # write.table(predict(cs.forest.grid, round(data$X, 2))$predictions, file = "data/causal_survival_predictions_grid.csv", row.names = FALSE, col.names = FALSE)

  expect_equal(predict(cs.forest.grid)$predictions, expected.predictions.oob.grid)
  expect_equal(predict(cs.forest.grid, round(data$X, 2))$predictions, expected.predictions.grid)
})

test_that("causal survival forest summary functions works as expected", {
  n <- 500
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, dgp = "simple1")
  X <- data$X
  W <- data$W
  Y <- data$Y
  D <- data$D
  tau <- data$cate
  cs.forest <- causal_survival_forest(X, Y, W, D, num.trees = 500)

  blp <- best_linear_projection(cs.forest)
  ate <- average_treatment_effect(cs.forest)
  blp.subset <- best_linear_projection(cs.forest, subset = X[, 1] > 0.5)
  ate.subset <- average_treatment_effect(cs.forest, subset = X[, 1] > 0.5)

  expect_equal(blp[1], ate[["estimate"]])
  expect_equal(blp[2], ate[["std.err"]], tolerance = 1e-4)
  expect_equal(ate[["estimate"]], mean(tau), tolerance = 3 * ate[["std.err"]])

  expect_equal(blp.subset[1], ate.subset[["estimate"]])
  expect_equal(blp.subset[2], ate.subset[["std.err"]], tolerance = 1e-4)
  expect_equal(ate.subset[["estimate"]], mean(tau[X[, 1] > 0.5]), tolerance = 3 * ate.subset[["std.err"]])
})
