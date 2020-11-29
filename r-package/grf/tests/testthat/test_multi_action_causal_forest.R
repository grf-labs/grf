library(grf)

test_that("single treatment multi_action_causal_forest is similar to causal_forest", {
  # It is not possible to check this parity holds exactly since forest differences
  # accrue through numerical differences (e.g. relabeling in causal forest is done with doubles
  # and with Eigen data structures in multi action causal forest.)
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  nmissing <- 50
  X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN

  cf <- causal_forest(X, Y, W, W.hat = 0, Y.hat = 0, seed = 1, stabilize.splits = FALSE,
                     alpha = 0, min.node.size = 1, num.trees = 500)
  mcf <- multi_action_causal_forest(X, Y, as.factor(W), W.hat = c(0, 0), Y.hat = 0, seed = 1,
                                    alpha = 0, min.node.size = 1, num.trees = 500)

  expect_equal(mean((predict(cf)$predictions - predict(mcf)$predictions[, 1])^2), 0, tol = 0.05)
  expect_equal(mean(predict(cf)$predictions), mean(predict(mcf)$predictions[, 1]), tol = 0.05)

  expect_equal(mean((predict(cf, X)$predictions - predict(mcf, X)$predictions[, 1])^2), 0, tol =  0.05)
  expect_equal(mean(predict(cf, X)$predictions), mean(predict(mcf, X)$predictions[, 1]), tol = 0.05)
})

test_that("multi_action_causal_forest contrasts works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  Y <- X[, 1] + 1.5 * (W == "A") + 2.8 * (W == "B") - 4 * (W == "C") + 0.1 * rnorm(n)

  mcf.A <- multi_action_causal_forest(X, Y, W, num.trees = 500, seed = 1)
  tau.hat.oob.A <- predict(mcf.A)$predictions
  tau.hat.A <- predict(mcf.A, X)$predictions

  mcf.C <- multi_action_causal_forest(X, Y, relevel(W, ref = "C"), num.trees = 500, seed = 1)
  tau.hat.oob.C <- predict(mcf.C)$predictions
  tau.hat.C <- predict(mcf.C, X)$predictions

  # 1. With easy constant treatment effects we estimate the correct contrasts
  expect_equal(colMeans(tau.hat.oob.A), c("B - A" = 2.8 - 1.5, "C - A" = -4 - 1.5), tol = 0.04)
  expect_equal(colMeans(tau.hat.A), c("B - A" = 2.8 - 1.5, "C - A" = -4 - 1.5), tol = 0.04)
  expect_equal(colMeans(tau.hat.oob.C), c("A - C" = 1.5 - (-4), "B - C" = 2.8 - (-4)), tol = 0.04)
  expect_equal(colMeans(tau.hat.C), c("A - C" = 1.5 - (-4), "B - C" = 2.8 - (-4)), tol = 0.04)

  # 2. The estimated contrast respects the symmetry properties we expect. It is not possible to check
  # this invariant exactly since differences in relabeling may lead to different trees
  expect_equal(tau.hat.oob.A[, "C - A"], -1 * tau.hat.oob.C[, "A - C"], tol = 0.05)
  expect_equal(tau.hat.A[, "C - A"], -1 * tau.hat.C[, "A - C"], tol = 0.05)

  expect_equal(tau.hat.oob.A[, "B - A"] - tau.hat.oob.A[, "C - A"], tau.hat.oob.C[, "B - C"], tol = 0.02)
  expect_equal(tau.hat.A[, "B - A"] - tau.hat.A[, "C - A"], tau.hat.C[, "B - C"], tol = 0.02)
})

test_that("multi_action_causal_forest ATE works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  tauB <- pmax(X[, 2], 0)
  tauC <- - 1.5 * abs(X[, 2])
  Y <- 2 + X[, 1] + tauB * (W == "B") + tauC * (W == "C") + 0.1 * rnorm(n)
  mcf <- multi_action_causal_forest(X, Y, W, num.trees = 500)

  ate <- average_treatment_effect(mcf)
  expect_equal(ate["B - A", "estimate"], mean(tauB), tol = 3 * ate["B - A", "std.err"])
  expect_equal(ate["C - A", "estimate"], mean(tauC), tol = 3 * ate["C - A", "std.err"])

  ate.subset <- average_treatment_effect(mcf, subset = X[, 2] < 0)
  expect_equal(ate.subset["B - A", "estimate"], 0, tol = 3 * ate.subset["B - A", "std.err"])

  mcf.B <- multi_action_causal_forest(X, Y, relevel(W, ref = "B"), num.trees = 500)
  ate.B <- average_treatment_effect(mcf.B)

  expect_equal(ate.B["A - B", "estimate"], -1 * ate["B - A", "estimate"], tol = 0.1)
  expect_equal(ate.B["C - B", "estimate"], ate["C - A", "estimate"] - ate["B - A", "estimate"], tol = 0.1)
})
