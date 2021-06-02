library(grf)

test_that("single treatment multi_arm_causal_forest is similar to causal_forest", {
  # It is not possible to check this parity holds exactly since forest differences
  # accrue through numerical differences (e.g. relabeling in causal forest is done with doubles
  # and with Eigen data structures in multi arm causal forest.)
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  nmissing <- 50
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  cf <- causal_forest(X, Y, W, W.hat = 1/2, Y.hat = 0, seed = 42, num.trees = 500)
  mcf <- multi_arm_causal_forest(X, Y, as.factor(W), W.hat = c(1/2, 1/2), Y.hat = 0, seed = 42, num.trees = 500)

  pp.cf <- predict(cf, estimate.variance = TRUE)
  pp.mcf <- predict(mcf, estimate.variance = TRUE)
  z.cf <- abs(pp.cf$predictions - tau) / sqrt(pp.cf$variance.estimates)
  z.mcf <- abs(pp.mcf$predictions[,,] - tau) / sqrt(pp.mcf$variance.estimates)
  expect_equal(mean(z.cf <= 1.96), mean(z.mcf <= 1.96), tolerance = 0.05)
  expect_equal(mean((pp.cf$predictions - pp.mcf$predictions)^2), 0, tolerance = 0.05)
  expect_equal(mean(pp.cf$predictions), mean(pp.mcf$predictions), tolerance = 0.05)
  expect_equal(mean((predict(cf, X)$predictions - predict(mcf, X)$predictions)^2), 0, tolerance = 0.05)
  expect_equal(mean(predict(cf, X)$predictions), mean(predict(mcf, X)$predictions), tolerance = 0.05)
  expect_equal(average_treatment_effect(cf)[["estimate"]], average_treatment_effect(mcf)[["estimate"]], tolerance = 0.001)
  expect_equal(average_treatment_effect(cf)[["std.err"]], average_treatment_effect(mcf)[["std.err"]], tolerance = 0.001)

  # Same checks with standard unconstrained regression splits.
  cf.rsplit <- causal_forest(X, Y, W, W.hat = 1/2, Y.hat = 0, seed = 42, num.trees = 500, stabilize.splits = FALSE)
  mcf.rsplit <- multi_arm_causal_forest(X, Y, as.factor(W), W.hat = c(1/2, 1/2), Y.hat = 0, seed = 42,
                                        num.trees = 500, stabilize.splits = FALSE)

  pp.cf.rsplit <- predict(cf.rsplit, estimate.variance = TRUE)
  pp.mcf.rsplit <- predict(mcf.rsplit, estimate.variance = TRUE)
  z.cf.rsplit <- abs(pp.cf.rsplit$predictions - tau) / sqrt(pp.cf.rsplit$variance.estimates)
  z.mcf.rsplit <- abs(pp.mcf.rsplit$predictions[,,] - tau) / sqrt(pp.mcf.rsplit$variance.estimates)
  expect_equal(mean(z.cf.rsplit <= 1.96), mean(z.mcf.rsplit <= 1.96), tolerance = 0.05)
  expect_equal(mean((pp.cf.rsplit$predictions - pp.mcf.rsplit$predictions)^2), 0, tolerance = 0.05)
  expect_equal(mean(pp.cf.rsplit$predictions), mean(pp.mcf.rsplit$predictions), tolerance = 0.05)
  expect_equal(mean((predict(cf.rsplit, X)$predictions - predict(mcf.rsplit, X)$predictions)^2), 0, tolerance = 0.05)
  expect_equal(mean(predict(cf.rsplit, X)$predictions), mean(predict(mcf.rsplit, X)$predictions), tolerance = 0.05)
  expect_equal(average_treatment_effect(cf.rsplit)[["estimate"]], average_treatment_effect(mcf.rsplit)[["estimate"]], tolerance = 0.001)
  expect_equal(average_treatment_effect(cf.rsplit)[["std.err"]], average_treatment_effect(mcf.rsplit)[["std.err"]], tolerance = 0.001)
})

test_that("multi_arm_causal_forest contrasts works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  Y <- X[, 1] + 1.5 * (W == "A") + 2.8 * (W == "B") - 4 * (W == "C") + 0.1 * rnorm(n)

  mcf.A <- multi_arm_causal_forest(X, Y, W, num.trees = 500, seed = 42)
  tau.hat.oob.A <- predict(mcf.A)$predictions[,,]
  tau.hat.A <- predict(mcf.A, X)$predictions[,,]

  mcf.C <- multi_arm_causal_forest(X, Y, relevel(W, ref = "C"), num.trees = 500, seed = 42)
  tau.hat.oob.C <- predict(mcf.C)$predictions[,,]
  tau.hat.C <- predict(mcf.C, X)$predictions[,,]

  # 1. With easy constant treatment effects we estimate the correct contrasts
  expect_equal(colMeans(tau.hat.oob.A), c("B - A" = 2.8 - 1.5, "C - A" = -4 - 1.5), tolerance = 0.04)
  expect_equal(colMeans(tau.hat.A), c("B - A" = 2.8 - 1.5, "C - A" = -4 - 1.5), tolerance = 0.04)
  expect_equal(colMeans(tau.hat.oob.C), c("A - C" = 1.5 - (-4), "B - C" = 2.8 - (-4)), tolerance = 0.04)
  expect_equal(colMeans(tau.hat.C), c("A - C" = 1.5 - (-4), "B - C" = 2.8 - (-4)), tolerance = 0.04)

  # 2. The estimated contrast respects the symmetry properties we expect. It is not possible to check
  # this invariant exactly since differences in relabeling may lead to different trees
  expect_equal(tau.hat.oob.A[, "C - A"], -1 * tau.hat.oob.C[, "A - C"], tolerance = 0.01)
  expect_equal(tau.hat.A[, "C - A"], -1 * tau.hat.C[, "A - C"], tolerance = 0.01)

  expect_equal(tau.hat.oob.A[, "B - A"] - tau.hat.oob.A[, "C - A"], tau.hat.oob.C[, "B - C"], tolerance = 0.01)
  expect_equal(tau.hat.A[, "B - A"] - tau.hat.A[, "C - A"], tau.hat.C[, "B - C"], tolerance = 0.01)

  # The above invariance holds exactly if we ignore splitting and just predict
  mcf.A.ns <- multi_arm_causal_forest(X, Y, W, num.trees = 250, seed = 42, min.node.size = n)
  tau.hat.oob.A.ns <- predict(mcf.A.ns)$predictions[,,]
  tau.hat.A.ns <- predict(mcf.A.ns, X)$predictions[,,]

  mcf.C.ns <- multi_arm_causal_forest(X, Y, relevel(W, ref = "C"), num.trees = 250, seed = 42, min.node.size = n)
  tau.hat.oob.C.ns <- predict(mcf.C.ns)$predictions[,,]
  tau.hat.C.ns <- predict(mcf.C.ns, X)$predictions[,,]

  expect_equal(tau.hat.oob.A.ns[, "C - A"], -1 * tau.hat.oob.C.ns[, "A - C"], tolerance = 1e-10)
  expect_equal(tau.hat.A.ns[, "C - A"], -1 * tau.hat.C.ns[, "A - C"], tolerance = 1e-10)

  expect_equal(tau.hat.oob.A.ns[, "B - A"] - tau.hat.oob.A.ns[, "C - A"], tau.hat.oob.C.ns[, "B - C"], tolerance = 1e-10)
  expect_equal(tau.hat.A.ns[, "B - A"] - tau.hat.A.ns[, "C - A"], tau.hat.C.ns[, "B - C"], tolerance = 1e-10)
})

test_that("multi_arm_causal_forest with binary treatment respects contrast invariance", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

  # With W.hat = 1/2 we can make the following check exact.
  # Setting W.hat to NULL or any other float pair will cause extremely small numerical differences
  # to accrue (10th+ digit) in splitting on the pseudo outcomes, leading these forests
  # to yield different splits even though they are algebraically the same.
  W.hat <- c(0.5, 0.5)
  cf <- multi_arm_causal_forest(X, Y, as.factor(W), W.hat = W.hat, num.trees = 250, seed = 42)
  cf.flipped <- multi_arm_causal_forest(X, Y, as.factor(1 - W), W.hat = W.hat, num.trees = 250, seed = 42)
  cf.relevel <- multi_arm_causal_forest(X, Y, relevel(as.factor(W), ref = "1"), W.hat = W.hat, num.trees = 250, seed = 42)

  pp <- predict(cf)$predictions[,,]
  pp.flipped <- predict(cf.flipped)$predictions[,,]
  pp.relevel <- predict(cf.relevel)$predictions[,,]

  expect_equal(pp, -1 * pp.flipped, tolerance = 0)
  expect_equal(pp, -1 * pp.relevel, tolerance = 0)
  expect_equal(pp.flipped, pp.relevel, tolerance = 0)
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
  expect_equal(ate["B - A", "estimate"], mean(tauB), tolerance = 3 * ate["B - A", "std.err"])
  expect_equal(ate["C - A", "estimate"], mean(tauC), tolerance = 3 * ate["C - A", "std.err"])

  ate.subset <- average_treatment_effect(mcf, subset = X[, 2] < 0)
  expect_equal(ate.subset["B - A", "estimate"], 0, tolerance = 3 * ate.subset["B - A", "std.err"])

  mcf.B <- multi_arm_causal_forest(X, Y, relevel(W, ref = "B"), num.trees = 500)
  ate.B <- average_treatment_effect(mcf.B)

  expect_equal(ate.B["A - B", "estimate"], -1 * ate["B - A", "estimate"], tolerance = 0.05)
  expect_equal(ate.B["C - B", "estimate"], ate["C - A", "estimate"] - ate["B - A", "estimate"], tolerance = 0.05)
})

test_that("multi_arm_causal_forest ATE standard errors are consistent with rest of GRF", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  tauB <- pmax(X[, 2], 0)
  tauC <- - 1.5 * abs(X[, 2])
  Y <- 2 + X[, 1] + tauB * (W == "B") + tauC * (W == "C") + rnorm(n)
  mcf <- multi_arm_causal_forest(X, Y, W, num.trees = 500)

  ate <- average_treatment_effect(mcf)
  DR.scores <- get_scores(mcf)

  lm.test <- lmtest::coeftest(lm(DR.scores[,,] ~ 1),
                              vcov = sandwich::vcovCL,
                              type = "HC3",
                              cluster = clusters <- if (length(mcf$clusters) > 0)
                                mcf$clusters else 1:length(mcf$Y.orig)
                              )
  expect_equal(ate[1, 2], lm.test[1, 2], tolerance = 1e-10)
  expect_equal(ate[2, 2], lm.test[2, 2], tolerance = 1e-10)
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
  alpha1 <- get_forest_weights(mcf, x1)[1, ]
  theta1.lm <- lm(Y ~ W.matrix[, -1], weights = alpha1)

  theta1.weighted <- predict(mcf.weighted, x1)$predictions
  alpha1.weighted <- get_forest_weights(mcf.weighted, x1)[1, ]
  theta1.lm.weighted <- lm(Y ~ W.matrix[, -1], weights = alpha1.weighted * sample.weights)

  expect_equal(as.numeric(theta1), as.numeric(theta1.lm$coefficients[-1]), tolerance = 1e-6)
  expect_equal(as.numeric(theta1.weighted), as.numeric(theta1.lm.weighted$coefficients[-1]), tolerance = 1e-6)
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
  forest.1 <- multi_arm_causal_forest(X, Y, W, sample.weights = sample.weights, num.trees = 250, seed = 42)
  forest.2 <- multi_arm_causal_forest(X, Y, W, sample.weights = 64 * sample.weights, num.trees = 250, seed = 42)
  pred.1 <- predict(forest.1, estimate.variance = TRUE)
  pred.2 <- predict(forest.2, estimate.variance = TRUE)

  expect_equal(pred.1$predictions, pred.2$predictions, tolerance = 1e-10)
  expect_equal(pred.1$variance.estimates, pred.2$variance.estimates, tolerance = 1e-10)
})

test_that("multi_arm_causal_forest confidence intervals are reasonable", {
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  tauB <- pmax(X[, 2], 0)
  tauC <- - 1.5 * abs(X[, 2])
  Y <- 2 + X[, 1] + tauB * (W == "B") + tauC * (W == "C") + rnorm(n)
  # Use min.node.size = 1 to check coverage on small n.
  mcf <- multi_arm_causal_forest(X, Y, W, num.trees = 500, min.node.size = 1)

  tau <- cbind(tauB, tauC)
  pp.mcf <- predict(mcf, estimate.variance = TRUE)
  z.score <- abs(pp.mcf$predictions[,,] - tau) / sqrt(pp.mcf$variance.estimates)
  coverage <- colMeans(z.score <= 1.96)

  expect_gt(coverage[1], 0.75)
  expect_gt(coverage[2], 0.75)
})

test_that("sample weighted multi_arm_causal forest gives correct coverage", {
  n <- 2500
  p <- 5
  pA <- 0.2
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, pA)
  gamma <- A / pA + (1 - A) / (1 - pA)
  tau.true <- 1 / (1 + exp(-3 * X[,1]))
  Y <- 2 * tau.true * W * A + rnorm(n)

  cf <- multi_arm_causal_forest(X, Y, as.factor(W), W.hat = c(0.5, 0.5), sample.weights = gamma)
  preds <- predict(cf, estimate.variance = TRUE)
  zstat <- (preds$predictions[,,] - tau.true) / sqrt(preds$variance.estimates)
  coverage <- mean(abs(zstat) < 1.96)
  mse <- mean((preds$predictions[,,] - tau.true)^2)

  cf.noweight <- multi_arm_causal_forest(X, Y, as.factor(W), W.hat = c(0.5, 0.5))
  preds.noweight <- predict(cf.noweight, estimate.variance = TRUE)
  zstat.noweight <- (preds.noweight$predictions[,,] - tau.true) / sqrt(preds.noweight$variance.estimates)
  coverage.noweight <- mean(abs(zstat.noweight) < 1.96)
  mse.noweight <- mean((preds.noweight$predictions[,,] - tau.true)^2)

  expect_lt(mse / mse.noweight, 0.3)
  expect_gt(coverage, 0.75)
  expect_lt(coverage.noweight, 0.55)
})

test_that("multi_arm_causal_forest with multiple outcomes works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  Y <- X[, 1] + X[, 2] * (W == "B") - 1.5 * X[, 2] * (W == "C") + rnorm(n)
  nmissing <- 50
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  # Check various symmetry properties.
  # A multi arm causal forest trained on Y is identical to one trained on [Y 0]
  rf <- multi_arm_causal_forest(X, Y, W, Y.hat = 0, W.hat = c(1/3, 1/3, 1/3),
                                num.trees = 200, seed = 42)
  ate.rf <- average_treatment_effect(rf)
  row.names(ate.rf) <- NULL
  mrf <- multi_arm_causal_forest(X, cbind(Y.1 = Y, 0), W, Y.hat = c(0, 0), W.hat = c(1/3, 1/3, 1/3),
                                 num.trees = 200, seed = 42)

  expect_equal(predict(rf)$predictions[,,], predict(mrf)$predictions[, , 1])
  expect_equal(predict(rf, X)$predictions[,,], predict(mrf, X)$predictions[, , 1])
  expect_equal(dim(predict(mrf)$predictions), c(n, 2, 2))
  ate.mrf <- average_treatment_effect(mrf)[1:2, ]
  row.names(mrf) <- NULL
  expect_equal(ate.rf, ate.mrf)

  # The above logic holds "symmetrically"
  # A multi arm causal forest trained on Y is identical to one trained on [0 Y]
  mrf <- multi_arm_causal_forest(X, cbind(0, Y.1 = Y), W, Y.hat = c(0, 0), W.hat = c(1/3, 1/3, 1/3),
                                 num.trees = 200, seed = 42)

  expect_equal(predict(rf)$predictions[,,], predict(mrf)$predictions[, , 2])
  expect_equal(predict(rf, X)$predictions[,,], predict(mrf, X)$predictions[, , 2])
  expect_equal(dim(predict(mrf)$predictions), c(n, 2, 2))
  ate.mrf <- average_treatment_effect(mrf)[3:4, ]
  row.names(ate.mrf) <- NULL
  expect_equal(ate.rf, ate.mrf)

  # A multi arm causal forest trained on Y is identical to one trained on [0 0 Y 0 0 0]
  mrf <- multi_arm_causal_forest(X, cbind(0, 0, Y.1 = Y, 0, 0, 0), W, Y.hat = rep(0, 6), W.hat = c(1/3, 1/3, 1/3),
                                 num.trees = 200, seed = 42)

  expect_equal(predict(rf)$predictions[,,], predict(mrf)$predictions[, , 3])
  expect_equal(predict(rf, X)$predictions[,,], predict(mrf, X)$predictions[, , 3])
  expect_equal(dim(predict(mrf)$predictions), c(n, 2, 6))
  ate.mrf <- average_treatment_effect(mrf)[5:6, ]
  row.names(ate.mrf) <- NULL
  expect_equal(ate.rf, ate.mrf)

  # A multi arm causal forest trained on duplicated Y's yields the same result
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  YY <- X[, 1] + X[, 2] * (W == "B") - 1.5 * X[, 2] * (W == "C") + matrix(rnorm(n * 2), n, 2)
  colnames(YY) <- 1:2
  mrf <- multi_arm_causal_forest(X, YY, W, Y.hat = c(0, 0), W.hat = c(1/3, 1/3, 1/3),
                                 num.trees = 200, seed = 42)
  mrf.dup <- multi_arm_causal_forest(X, cbind(YY, YY), W, Y.hat = rep(0, 4), W.hat = c(1/3, 1/3, 1/3),
                                     num.trees = 200, seed = 42)

  expect_equal(unname(predict(mrf)$predictions), unname(predict(mrf.dup)$predictions[, , 1:2]))
  expect_equal(unname(predict(mrf)$predictions), unname(predict(mrf.dup)$predictions[, , 3:4]))
  expect_equal(average_treatment_effect(mrf), average_treatment_effect(mrf.dup)[1:4, ])
  df1 <- average_treatment_effect(mrf)[, -4]
  df2 <- average_treatment_effect(mrf.dup)[5:8, -4]
  row.names(df2) <- NULL
  expect_equal(df1, df2)
})

test_that("multi_arm_causal_forest with multiple outcomes is well calibrated", {
  # Simulate n correlated mean zero normal draws with covariance matrix sigma
  # using the Cholesky decomposition. Returns a [n X ncol(sigma)] matrix.
  rmvnorm <- function(n, sigma) {
    K <- ncol(sigma)
    A <- chol(sigma)
    z <- matrix(rnorm(n * K), n, K)
    z %*% A
  }

  # A multi arm causal forest fit on two outcomes yields lower MSE than two separate
  # causal forests when the CATEs are correlated with low idiosyncratic noise.
  sigma <- diag(4)
  sigma[2, 1] <- sigma[1, 2] <- 0.5
  sigma[3, 4] <- sigma[4, 3] <- 0.5

  n <- 500
  p <- 4
  X <- rmvnorm(n, sigma = sigma)
  W <- rbinom(n, 1, 0.5)
  tau1 <- pmax(X[, 1], 0)
  tau2 <- pmax(X[, 2], 0)
  tau <- cbind(tau1, tau2)
  YY <- tau * W + X[, 3:4] + 0.5 * matrix(rnorm(n * 2), n, 2)

  cf.pred <- apply(YY, 2, function(Y) predict(causal_forest(X, Y, W, num.trees = 500))$predictions)
  mcf.pred <- predict(multi_arm_causal_forest(X, YY, as.factor(W), num.trees = 500))$predictions[,,]
  mse.cf <- mean((mcf.pred - tau)^2)
  mse.mcf <- mean((cf.pred - tau)^2)

  expect_lt(mse.mcf / mse.cf,  0.85)
})
