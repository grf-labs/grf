test_that("causal survival forest is well-calibrated", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, n.mc = 5000, dgp = "simple1")
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 500)
  cs.forest.prob <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 500,
                                           target = "survival.probability", horizon = data$y0)
  cs.pred <- predict(cs.forest)
  cs.pred.prob <- predict(cs.forest.prob)
  mse.oob <- mean((cs.pred$predictions - data$cate)^2)
  mse.oob.prob <- mean((cs.pred.prob$predictions - data$cate.prob)^2)

  X.test <- matrix(0.5, 10, p)
  X.test[, 1] <- seq(0, 1, length.out = 10)
  cs.pred.test <- predict(cs.forest, X.test)
  cs.pred.test.prob <- predict(cs.forest.prob, X.test)
  true.effect.test <- generate_causal_survival_data(10, p, Y.max = 1, X = X.test, n.mc = 5000, dgp = "simple1")
  mse.test <- mean((cs.pred.test$predictions - true.effect.test$cate)^2)
  mse.prob.test <- mean((cs.pred.test.prob$predictions - true.effect.test$cate.prob)^2)

  expect_lt(mse.oob, 0.01)
  expect_lt(mse.test, 0.01)
  expect_lt(mse.oob.prob, 0.02)
  expect_lt(mse.prob.test, 0.02)
})

test_that("causal survival forest with complete-data is ~identical to causal forest", {
  n <- 500
  p <- 5
  X <- matrix(runif(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y.max <- 1
  failure.time <- pmin(rexp(n) * X[, 1] + W, Y.max)
  censor.time <- 9999
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)

  # RMST
  cs.forest <- causal_survival_forest(X, Y, W, D, horizon = Y.max, num.trees = 500)
  pp.cs <- predict(cs.forest, estimate.variance = TRUE)
  cf <- causal_forest(X, Y, W, num.trees = 500)
  pp.cf <- predict(cf, estimate.variance = TRUE)
  ate.cs <- average_treatment_effect(cs.forest)
  ate.cf <- average_treatment_effect(cf)
  muhat0 <- cf$Y.hat - cf$W.hat * pp.cf$predictions
  muhat1 <- cf$Y.hat + (1 - cf$W.hat) * pp.cf$predictions
  muhat0.cs <- cs.forest$Y.hat - cs.forest$W.hat * pp.cs$predictions
  muhat1.cs <- cs.forest$Y.hat + (1 - cs.forest$W.hat) * pp.cs$predictions

  expect_lt(mean((pp.cs$predictions - pp.cf$predictions)^2), 0.001)
  expect_equal(mean(pp.cs$predictions), mean(pp.cf$predictions), tolerance = 0.02)
  expect_lt(mean((pp.cs$variance.estimates - pp.cf$variance.estimates)^2), 0.001)
  expect_equal(mean(pp.cs$variance.estimates), mean(pp.cf$variance.estimates), tolerance = 0.01)
  expect_equal(ate.cs[["estimate"]], ate.cf[["estimate"]], tolerance = 0.01)
  expect_equal(ate.cs[["std.err"]], ate.cf[["std.err"]], tolerance = 0.005)
  expect_equal(muhat0, muhat0.cs, tolerance = 0.15)
  expect_equal(muhat1, muhat1.cs, tolerance = 0.15)

  # Probability
  horizon <- 0.5
  cs.forest.prob <- causal_survival_forest(X, Y, W, D, target = "survival.probability", horizon = horizon, num.trees = 500)
  pp.cs.prob <- predict(cs.forest.prob, estimate.variance = TRUE)
  cf.prob <- causal_forest(X, as.numeric(Y > horizon), W, num.trees = 500)
  pp.cf.prob <- predict(cf.prob, estimate.variance = TRUE)
  ate.cs.prob <- average_treatment_effect(cs.forest.prob)
  ate.cf.prob <- average_treatment_effect(cf.prob)
  muhat0.prob <- cf.prob$Y.hat - cf.prob$W.hat * pp.cf.prob$predictions
  muhat1.prob <- cf.prob$Y.hat + (1 - cf.prob$W.hat) * pp.cf.prob$predictions
  muhat0.cs.prob <- cs.forest.prob$Y.hat - cs.forest.prob$W.hat * pp.cs.prob$predictions
  muhat1.cs.prob <- cs.forest.prob$Y.hat + (1 - cs.forest.prob$W.hat) * pp.cs.prob$predictions

  expect_lt(mean((pp.cs.prob$predictions - pp.cf.prob$predictions)^2), 0.001)
  expect_equal(mean(pp.cs.prob$predictions), mean(pp.cf.prob$predictions), tolerance = 0.02)
  expect_lt(mean((pp.cs.prob$variance.estimates - pp.cf.prob$variance.estimates)^2), 0.001)
  expect_equal(mean(pp.cs.prob$variance.estimates), mean(pp.cf.prob$variance.estimates), tolerance = 0.01)
  expect_equal(ate.cs.prob[["estimate"]], ate.cf.prob[["estimate"]], tolerance = 0.01)
  expect_equal(ate.cs.prob[["std.err"]], ate.cf.prob[["std.err"]], tolerance = 0.005)
  expect_equal(muhat0.prob, muhat0.cs.prob, tolerance = 0.15)
  expect_equal(muhat1.prob, muhat1.cs.prob, tolerance = 0.15)
})

test_that("causal survival forest probability predictions are close for close horizon", {
  n <- 500
  p <- 5
  data <- generate_causal_survival_data(n, p, dgp = "simple1")
  X <- data$X
  W <- data$W
  Y <- round(data$Y, 2)
  D <- data$D

  horizon <- median(Y)
  cs.forest <- causal_survival_forest(X, Y, W, D, target = "RMST", horizon = horizon, num.trees = 500)
  cs.forest.shiftp <- causal_survival_forest(X, Y, W, D, target = "RMST", horizon = horizon + 0.01, num.trees = 500)
  cs.forest.shiftn <- causal_survival_forest(X, Y, W, D, target = "RMST", horizon = horizon - 0.01, num.trees = 500)
  expect_equal(predict(cs.forest)$predictions, predict(cs.forest.shiftp)$predictions, tolerance = 0.05)
  expect_equal(predict(cs.forest)$predictions, predict(cs.forest.shiftn)$predictions, tolerance = 0.05)
})

test_that("causal survival forest is invariant to shifting all event times", {
  n <- 500
  p <- 5
  data <- generate_causal_survival_data(n, p, dgp = "simple1")
  X <- data$X
  W <- data$W
  Y <- round(data$Y, 2)
  D <- data$D

  cs.forest <- causal_survival_forest(X, Y, W, D, horizon = data$Y.max, seed = 42, num.trees = 500)
  cs.forest.shift <- causal_survival_forest(X, Y + 10, W, D, horizon = max(Y + 10), seed = 42, num.trees = 500)
  expect_equal(predict(cs.forest)$predictions, predict(cs.forest.shift)$predictions, tolerance = 1e-14)

  cs.forest.prob <- causal_survival_forest(X, Y, W, D, target = "survival.probability", horizon = 0.5, seed = 42, num.trees = 500)
  cs.forest.shift.prob <- causal_survival_forest(X, Y + 10, W, D, target = "survival.probability", horizon = 0.5 + 10, seed = 42, num.trees = 500)
  expect_equal(predict(cs.forest.prob)$predictions, predict(cs.forest.shift.prob)$predictions, tolerance = 1e-14)
})

test_that("causal survival forest predictions are kernel weighted correctly", {
  # Test that the sufficient statistic approach to solving the forest weighted
  # estimating equation is internally consistent.
  n <- 250
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, n.mc = 1, dgp = "simple1")
  sample.weights <- sample(c(1, 10), n, TRUE)
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 250)
  cs.forest.weighted <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max,
                                               num.trees = 250, sample.weights = sample.weights)
  x1 <- data$X[1, , drop = FALSE]
  cs.pred <- predict(cs.forest, x1)$predictions
  cs.pred.weighted <- predict(cs.forest.weighted, x1)$predictions
  forest.weights <- get_forest_weights(cs.forest, x1)[1, ]
  forest.weights.weighted <- get_forest_weights(cs.forest.weighted, x1)[1, ]

  theta1 <- sum(forest.weights * cs.forest[["_psi"]]$numerator) /
    sum(forest.weights * cs.forest[["_psi"]]$denominator)
  theta1.weighted <- sum(forest.weights.weighted * sample.weights * cs.forest.weighted[["_psi"]]$numerator) /
    sum(forest.weights.weighted * sample.weights * cs.forest.weighted[["_psi"]]$denominator)

  expect_equal(cs.pred, theta1)
  expect_equal(cs.pred.weighted, theta1.weighted)
})

test_that("causal survival forest variance estimates are decent", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, n.mc = 5000, dgp = "simple1")
  true.effect <- data$cate
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 500)
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
  cs.forest.weighted <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max,
                                               num.trees = 500, sample.weights = sample.weights)
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
  sf <- causal_survival_forest(X, Y, W, D, horizon = max(Y), num.trees = 100, sample.weights = weights, seed = 1)
  sf2 <- causal_survival_forest(X, Y, W, D, horizon = max(Y), num.trees = 100, sample.weights = weights * 2, seed = 1)
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
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 500)
  cs.pred <- predict(cs.forest)

  Y.grid <- seq(min(data$Y), max(data$Y), length.out = 25)
  cs.forest.grid <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 500,
                                           failure.times = Y.grid)
  cs.pred.grid <- predict(cs.forest.grid)

  expect_lt(mean((cs.pred$predictions - cs.pred.grid$predictions)^2), 0.001)
})

test_that("causal survival forest works as expected with missing values", {
  n <- 1000
  p <- 5
  data <- generate_causal_survival_data(n, p, n.mc = 1, dgp = "simple1")
  nmissing <- 200
  data$X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  csf <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 500)

  # MIA with data duplication
  Xl <- data$X
  Xr <- data$X
  Xl[is.nan(Xl)] <- -1e9
  Xr[is.nan(Xr)] <- 1e9
  X.mia <- cbind(Xl, Xr)
  csf.mia <- causal_survival_forest(X.mia, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 500)

  mse.oob.diff <- mean((predict(csf)$predictions - predict(csf.mia)$predictions)^2)
  mse.diff <- mean((predict(csf, data$X)$predictions - predict(csf.mia, X.mia)$predictions)^2)

  expect_equal(mse.oob.diff, 0, tolerance = 0.001)
  expect_equal(mse.diff, 0, tolerance = 0.001)
})

test_that("causal survival forest with survival target works as expected", {
  n <- 500
  p <- 5
  data <- generate_causal_survival_data(n, p, n.mc = 1, dgp = "simple1")
  cs.default <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 250,
                                       target = "survival.probability", horizon = 0.4, seed = 42)
  cs.full.grid <- causal_survival_forest(data$X, data$Y, data$W, data$D, num.trees = 250,
                                         failure.times = sort(unique(data$Y)),
                                         target = "survival.probability", horizon = 0.4, seed = 42)
  expect_equal(predict(cs.full.grid)$predictions, predict(cs.default)$predictions)
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
  tau.prob <- data$cate.prob
  cs.forest <- causal_survival_forest(X, Y, W, D, horizon = data$Y.max, num.trees = 500)
  cs.forest.prob <- causal_survival_forest(X, Y, W, D, target = "survival.probability",
                                           horizon = data$y0, num.trees = 500)

  blp <- best_linear_projection(cs.forest)
  ate <- average_treatment_effect(cs.forest)
  ate.prob <- average_treatment_effect(cs.forest.prob)
  blp.subset <- best_linear_projection(cs.forest, subset = X[, 1] > 0.5)
  ate.subset <- average_treatment_effect(cs.forest, subset = X[, 1] > 0.5)

  expect_equal(blp[1], ate[["estimate"]])
  expect_equal(blp[2], ate[["std.err"]], tolerance = 1e-4)
  expect_equal(ate[["estimate"]], mean(tau), tolerance = 3 * ate[["std.err"]])
  expect_equal(ate.prob[["estimate"]], mean(tau.prob), tolerance = 3 * ate.prob[["std.err"]])

  expect_equal(blp.subset[1], ate.subset[["estimate"]])
  expect_equal(blp.subset[2], ate.subset[["std.err"]], tolerance = 1e-4)
  expect_equal(ate.subset[["estimate"]], mean(tau[X[, 1] > 0.5]), tolerance = 3 * ate.subset[["std.err"]])
})

test_that("causal survival forest for partial effect estimation works as expected", {
  n <- 500
  p <- 5
  data <- generate_causal_survival_data(n, p, Y.max = 1, dgp = "simple1")
  X <- data$X
  W <- data$W
  Y <- data$Y
  D <- data$D
  cs.forest <- causal_survival_forest(X, Y, W, D, horizon = data$Y.max, num.trees = 500)
  cs.forest.prob <- causal_survival_forest(X, Y, W, D, target = "survival.probability",
                                           horizon = data$y0, num.trees = 500)
  # The following forests should be ~ identical
  W[1] <- W[1] + 1e-10
  cs.forest.ape <- causal_survival_forest(X, Y, W, D, horizon = data$Y.max, num.trees = 500)
  cs.forest.prob.ape <- causal_survival_forest(X, Y, W, D, target = "survival.probability",
                                               horizon = data$y0, num.trees = 500)
  expect_equal(cs.forest$predictions, cs.forest.ape$predictions, tolerance = 0.035)
  expect_equal(cs.forest.prob$predictions, cs.forest.prob.ape$predictions, tolerance = 0.035)

  ate <- average_treatment_effect(cs.forest)
  ate.ape <- average_treatment_effect(cs.forest.ape)
  ate.prob <- average_treatment_effect(cs.forest.prob)
  ate.prob.ape <- average_treatment_effect(cs.forest.prob.ape)
  expect_equal(ate[1], ate.ape[1], tolerance = 0.015)
  expect_equal(ate[2], ate.ape[2], tolerance = 0.0015)
  expect_equal(ate.prob[1], ate.prob.ape[1], tolerance = 0.015)
  expect_equal(ate.prob[2], ate.prob.ape[2], tolerance = 0.0015)

  blp <- best_linear_projection(cs.forest)
  blp.ape <- best_linear_projection(cs.forest.ape)
  blp.prob <- best_linear_projection(cs.forest.prob)
  blp.prob.ape <- best_linear_projection(cs.forest.prob.ape)
  expect_equal(blp[1], blp.ape[1], tolerance = 0.015)
  expect_equal(blp[2], blp.ape[2], tolerance = 0.0015)
  expect_equal(blp.prob[1], blp.prob.ape[1], tolerance = 0.015)
  expect_equal(blp.prob[2], blp.prob.ape[2], tolerance = 0.0015)

  # CSF estimates what it's supposed to in a simple setup
  n <- 1500
  p <- 5
  X <- matrix(runif(n*p), n, p)
  W <- runif(n, max = X[, 2] * 2)
  tau <- pmax(X[, 1], 0.5)
  Yt <- tau * W + X[, 4] + runif(n)
  censor.time <- 5 * runif(n)
  Y <- round(pmin(Yt, censor.time), 1)
  D <- as.integer(Yt <= censor.time)

  csf <- suppressWarnings(causal_survival_forest(X, Y, W, D, horizon = max(Y)))
  ate <- average_treatment_effect(csf)
  expect_equal(ate[["estimate"]], mean(tau), tolerance = 3 * ate[["std.err"]])

  # a complete-data CSF is ~same as CF
  p <- 6
  n <- 1000
  X <- matrix(2 * runif(n * p) - 1, n, p)
  W <- rbinom(n, 1, 0.5) + abs(rnorm(n))
  TAU <- 4 * (X[, 1] > 0)
  Y <- TAU * W + runif(n)

  cf <- causal_forest(X, Y, W, num.trees = 500)
  csf.complete <- causal_survival_forest(X, Y, W, rep(1, n), horizon = max(Y), num.trees = 500)
  ate.cf <- average_treatment_effect(cf)
  ate.csf.complete <- average_treatment_effect(csf.complete)
  expect_equal(ate.cf[[1]], ate.csf.complete[[1]], tolerance = 0.02)
  expect_equal(ate.cf[[2]], ate.csf.complete[[2]], tolerance = 0.005)
})

test_that("causal survival forest utility functions are internally consistent", {
  n <- 500
  p <- 5
  X <- matrix(runif(n * p), n, p)
  failure.time <- 0.25
  censor.time <- 0.5 * runif(n)
  Y <- pmin(failure.time, censor.time)
  D <- as.integer(failure.time <= censor.time)

  sf <- survival_forest(X, Y, D, num.trees = 500)
  pp.sf <- predict(sf)
  Y.hat <- expected_survival(pp.sf$predictions, pp.sf$failure.times)

  expect_equal(Y.hat, rep(failure.time, n))
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
  cs.forest <- causal_survival_forest(round(data$X, 2), round(data$Y, 2), data$W, data$D, horizon = data$Y.max,
                                      num.trees = 50, seed = 42, num.threads = 4)

  # Update with:
  # write.table(predict(cs.forest)$predictions, file = "data/causal_survival_oob_predictions.csv", row.names = FALSE, col.names = FALSE)
  # write.table(predict(cs.forest, round(data$X, 2))$predictions, file = "data/causal_survival_predictions.csv", row.names = FALSE, col.names = FALSE)
  expected.predictions.oob <- as.numeric(readLines("data/causal_survival_oob_predictions.csv"))
  expected.predictions <- as.numeric(readLines("data/causal_survival_predictions.csv"))

  expect_equal(predict(cs.forest)$predictions, expected.predictions.oob)
  expect_equal(predict(cs.forest, round(data$X, 2))$predictions, expected.predictions)

  # With target = "survival.probability"
  cs.forest.prob <- causal_survival_forest(round(data$X, 2), round(data$Y, 2), data$W, data$D,
                                           target = "survival.probability", horizon = 0.5,
                                           num.trees = 50, seed = 42, num.threads = 4)

  # Update with:
  # write.table(predict(cs.forest.prob)$predictions, file = "data/causal_survival_oob_predictions_prob.csv", row.names = FALSE, col.names = FALSE)
  # write.table(predict(cs.forest.prob, round(data$X, 2))$predictions, file = "data/causal_survival_predictions_prob.csv", row.names = FALSE, col.names = FALSE)
  expected.predictions.oob.prob <- as.numeric(readLines("data/causal_survival_oob_predictions_prob.csv"))
  expected.predictions.prob <- as.numeric(readLines("data/causal_survival_predictions_prob.csv"))

  expect_equal(predict(cs.forest.prob)$predictions, expected.predictions.oob.prob)
  expect_equal(predict(cs.forest.prob, round(data$X, 2))$predictions, expected.predictions.prob)
})
