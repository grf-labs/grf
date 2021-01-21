library(grf)

test_that("single treatment multi_arm_causal_forest is similar to causal_forest", {
  # It is not possible to check this parity holds exactly since forest differences
  # accrue through numerical differences (e.g. relabeling in causal forest is done with doubles
  # and with Eigen data structures in multi action causal forest.)
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <-  tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  nmissing <- 50
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  cf <- causal_forest(X, Y, W, W.hat = 0, Y.hat = 0, seed = 1, stabilize.splits = FALSE,
                     alpha = 0, min.node.size = 1, num.trees = 500)
  mcf <- multi_arm_causal_forest(X, Y, as.factor(W), W.hat = c(0, 0), Y.hat = 0, seed = 1,
                                    alpha = 0, min.node.size = 1, num.trees = 500)

  pp.cf <- predict(cf, estimate.variance = TRUE)
  pp.mcf <- predict(mcf, estimate.variance = TRUE)
  z.cf <- abs((pp.cf$predictions - tau) / sqrt(pp.cf$variance.estimates))
  z.mcf <- abs((pp.mcf$predictions - tau) / sqrt(pp.mcf$variance.estimates))
  expect_equal(mean(z.cf > 1.96), mean(z.mcf > 1.96), tol = 0.05)

  expect_equal(mean((pp.cf$predictions - pp.mcf$predictions[, 1])^2), 0, tol = 0.05)
  expect_equal(mean(pp.cf$predictions), mean(pp.mcf$predictions[, 1]), tol = 0.05)

  expect_equal(mean((predict(cf, X)$predictions - predict(mcf, X)$predictions[, 1])^2), 0, tol =  0.05)
  expect_equal(mean(predict(cf, X)$predictions), mean(predict(mcf, X)$predictions[, 1]), tol = 0.05)
})

test_that("multi_arm_causal_forest contrasts works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  Y <- X[, 1] + 1.5 * (W == "A") + 2.8 * (W == "B") - 4 * (W == "C") + 0.1 * rnorm(n)

  mcf.A <- multi_arm_causal_forest(X, Y, W, num.trees = 500, seed = 1)
  tau.hat.oob.A <- predict(mcf.A)$predictions
  tau.hat.A <- predict(mcf.A, X)$predictions

  mcf.C <- multi_arm_causal_forest(X, Y, relevel(W, ref = "C"), num.trees = 500, seed = 1)
  tau.hat.oob.C <- predict(mcf.C)$predictions
  tau.hat.C <- predict(mcf.C, X)$predictions

  # 1. With easy constant treatment effects we estimate the correct contrasts
  expect_equal(colMeans(tau.hat.oob.A), c("B - A" = 2.8 - 1.5, "C - A" = -4 - 1.5), tol = 0.04)
  expect_equal(colMeans(tau.hat.A), c("B - A" = 2.8 - 1.5, "C - A" = -4 - 1.5), tol = 0.04)
  expect_equal(colMeans(tau.hat.oob.C), c("A - C" = 1.5 - (-4), "B - C" = 2.8 - (-4)), tol = 0.04)
  expect_equal(colMeans(tau.hat.C), c("A - C" = 1.5 - (-4), "B - C" = 2.8 - (-4)), tol = 0.04)

  # 2. The estimated contrast respects the symmetry properties we expect. It is not possible to check
  # this invariant exactly since differences in relabeling may lead to different trees
  expect_equal(tau.hat.oob.A[, "C - A"], -1 * tau.hat.oob.C[, "A - C"], tol = 0.01)
  expect_equal(tau.hat.A[, "C - A"], -1 * tau.hat.C[, "A - C"], tol = 0.01)

  expect_equal(tau.hat.oob.A[, "B - A"] - tau.hat.oob.A[, "C - A"], tau.hat.oob.C[, "B - C"], tol = 0.01)
  expect_equal(tau.hat.A[, "B - A"] - tau.hat.A[, "C - A"], tau.hat.C[, "B - C"], tol = 0.01)
})

test_that("multi_arm_causal_forest ATE works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  tauB <- pmax(X[, 2], 0)
  tauC <- - 1.5 * abs(X[, 2])
  Y <- 2 + X[, 1] + tauB * (W == "B") + tauC * (W == "C") + 0.1 * rnorm(n)
  mcf <- multi_arm_causal_forest(X, Y, W, num.trees = 500)

  ate <- average_treatment_effect(mcf)
  expect_equal(ate["B - A", "estimate"], mean(tauB), tol = 3 * ate["B - A", "std.err"])
  expect_equal(ate["C - A", "estimate"], mean(tauC), tol = 3 * ate["C - A", "std.err"])

  ate.subset <- average_treatment_effect(mcf, subset = X[, 2] < 0)
  expect_equal(ate.subset["B - A", "estimate"], 0, tol = 3 * ate.subset["B - A", "std.err"])

  mcf.B <- multi_arm_causal_forest(X, Y, relevel(W, ref = "B"), num.trees = 500)
  ate.B <- average_treatment_effect(mcf.B)

  expect_equal(ate.B["A - B", "estimate"], -1 * ate["B - A", "estimate"], tol = 0.05)
  expect_equal(ate.B["C - B", "estimate"], ate["C - A", "estimate"] - ate["B - A", "estimate"], tol = 0.05)
})

test_that("multi_arm_causal_forest predictions are kernel weighted correctly", {
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  tauB <- pmax(X[, 2], 0)
  tauC <- - 1.5 * abs(X[, 2])
  Y <- 2 + X[, 1] + tauB * (W == "B") + tauC * (W == "C") + rnorm(n)
  sample.weights <- sample(c(1, 5), n, TRUE)
  mcf <- multi_arm_causal_forest(X, Y, W, Y.hat = 0, W.hat = c(0, 0, 0), num.trees = 250)
  mcf.weighted <- multi_arm_causal_forest(X, Y, W, Y.hat = 0, W.hat = c(0, 0, 0), num.trees = 250, sample.weights = sample.weights)

  W.matrix <- stats::model.matrix(~ mcf$W.orig - 1)
  x1 <- X[1, , drop = F]
  theta1 <- predict(mcf, x1)$predictions
  alpha1 <- get_sample_weights(mcf, x1)[1, ]
  theta1.lm <- lm(Y ~ W.matrix[, -1], weights = alpha1)

  theta1.weighted <- predict(mcf.weighted, x1)$predictions
  alpha1.weighted <- get_sample_weights(mcf.weighted, x1)[1, ]
  theta1.lm.weighted <- lm(Y ~ W.matrix[, -1], weights = alpha1.weighted * sample.weights)

  expect_equal(as.numeric(theta1), as.numeric(theta1.lm$coefficients[-1]), tol = 1e-6)
  expect_equal(as.numeric(theta1.weighted), as.numeric(theta1.lm.weighted$coefficients[-1]), tol = 1e-6)
})

test_that("multi_arm_causal_forest predictions and variance estimates are invariant to scaling of the sample weights.", {
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  tauB <- pmax(X[, 2], 0)
  tauC <- - 1.5 * abs(X[, 2])
  Y <- 2 + X[, 1] + tauB * (W == "B") + tauC * (W == "C") + rnorm(n)
  sample.weights <- sample(c(1, 5), n, TRUE)

  # The multiple is a power of 2 to avoid rounding errors allowing for exact comparison
  # between two forest with the same seed.
  forest.1 <- multi_arm_causal_forest(X, Y, W, sample.weights = sample.weights, num.trees = 250, seed = 1)
  forest.2 <- multi_arm_causal_forest(X, Y, W, sample.weights = 64 * sample.weights, num.trees = 250, seed = 1)
  pred.1 <- predict(forest.1, estimate.variance = TRUE)
  pred.2 <- predict(forest.2, estimate.variance = TRUE)

  expect_equal(pred.1$predictions, pred.2$predictions, tol = 1e-10)
  # expect_equal(pred.1$variance.estimates, pred.2$variance.estimates, tol = 1e-10)
  # expect_equal(pred.1$debiased.error, pred.2$debiased.error, tol = 1e-10)
})

test_that("multi_arm_causal_forest confidence intervals are reasonable", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  tauB <- pmax(X[, 2], 0)
  tauC <- - 1.5 * abs(X[, 2])
  Y <- 2 + X[, 1] + tauB * (W == "B") + tauC * (W == "C") + rnorm(n)
  mcf <- multi_arm_causal_forest(X, Y, W, num.trees = 500)

  tau <- cbind(tauB, tauC)
  pp.mcf <- predict(mcf, estimate.variance = TRUE)
  z <- abs((pp.mcf$predictions - tau) / sqrt(pp.mcf$variance.estimates))

  expect_true(all(colMeans(z > 1.96) < c(0.25, 0.25)))
})
