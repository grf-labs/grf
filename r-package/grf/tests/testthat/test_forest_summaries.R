library(grf)

set.seed(1234)

test_that("causal forest calibration is reasonable", {
  n <- 800
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 1], 0) * (W - 0.75) + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.25 + 0.5 * (X[, 1] > 0),
    Y.hat = pmax(X[, 1], 0) * (0.5 * (X[, 1] > 0) - 0.5),
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 1] - 1) <= 0.4)
  expect_true(abs(tc[2, 1] - 1) <= 0.4)
})

test_that("causal forest calibration is reasonable with no average effect", {
  n <- 800
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rnorm(n, 1, 0.5)
  Y <- sign(X[, 1]) * (W - 0.5) + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.5,
    Y.hat = 0,
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 3]) <= 3)
  expect_true(abs(tc[2, 1] - 1) <= 0.3)
})

test_that("causal forest calibration is reasonable with no heterogeneous effect", {
  n <- 800
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 2], 0) + W + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.25 + 0.5 * (X[, 1] > 0),
    Y.hat = 0.25 + 0.5 * (X[, 1] > 0) + pmax(X[, 2], 0),
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 1] - 1) <= 0.3)
  expect_true(abs(tc[2, 3]) <= 4)
})

test_that("causal forest calibration is reasonable with no heterogeneous effect with sample weights and clusters", {
  p <- 4
  K <- 100
  cluster.sizes <- pmax(20, round(40 + 3 * rt(K, df = 3)))
  n <- sum(cluster.sizes)
  clust <- rep(1:K, cluster.sizes)

  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  big <- as.numeric(cluster.sizes[clust] >= median(cluster.sizes))
  Y <- pmax(X[, 2], 0) * big + W + rnorm(n)

  e.cc <- 1 / (1 + exp(-2 * X[, 1]))
  cc <- as.logical(rbinom(n, 1, e.cc))
  ipcc.sample.weights <- 1 / e.cc
  clust.sample.weights <- 1 / cluster.sizes[clust]
  sample.weights <- ipcc.sample.weights * clust.sample.weights

  cf <- causal_forest(X[cc, ], Y[cc], W[cc],
    W.hat = 0.25 + 0.5 * (X[cc, 1] > 0),
    Y.hat = 0.25 + 0.5 * (X[cc, 1] > 0) + pmax(X[cc, 2], 0) * big[cc],
    sample.weights = sample.weights[cc],
    clusters = clust[cc],
    equalize.cluster.weights = FALSE,
    num.trees = 500
  )
  tc <- test_calibration(cf)

  expect_true(abs(tc[1, 1] - 1) <= 0.3)
  expect_true(abs(tc[2, 3]) <= 4)
})

test_that("regression forest calibration is reasonable", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- 5 + 5 * sign(X[, 1]) + rnorm(n)

  rf <- regression_forest(X, Y)
  tc <- test_calibration(rf)

  expect_true(abs(tc[1, 1] - 1) <= 0.1)
  expect_true(abs(tc[2, 1] - 1) <= 0.2)
})

test_that("regression forest calibration is reasonable with no heterogeneous effect", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- 5 + rnorm(n)

  rf <- regression_forest(X, Y)
  tc <- test_calibration(rf)

  expect_true(abs(tc[1, 1] - 1) <= 0.1)
  expect_true(abs(tc[2, 3]) <= 4)
})

test_that("causal forest calibration works with clusters", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 2], 0) + W + rnorm(n)

  cf <- causal_forest(X, Y, W,
    W.hat = 0.25 + 0.5 * (X[, 1] > 0),
    Y.hat = 0.25 + 0.5 * (X[, 1] > 0) + pmax(X[, 2], 0),
    num.trees = 100
  )
  tc <- test_calibration(cf)

  cf.clust <- cf
  cf.clust$W.orig <- c(cf$W.orig[1:(n / 2)], rep(cf$W.orig[n / 2 + 1:(n / 2)], 10))
  cf.clust$Y.orig <- c(cf$Y.orig[1:(n / 2)], rep(cf$Y.orig[n / 2 + 1:(n / 2)], 10))
  cf.clust$W.hat <- c(cf$W.hat[1:(n / 2)], rep(cf$W.hat[n / 2 + 1:(n / 2)], 10))
  cf.clust$Y.hat <- c(cf$Y.hat[1:(n / 2)], rep(cf$Y.hat[n / 2 + 1:(n / 2)], 10))
  cf.clust$predictions <- c(cf$predictions[1:(n / 2)], rep(cf$predictions[n / 2 + 1:(n / 2)], 10))
  cf.clust$debiased.error <- c(cf$debiased.error[1:(n / 2)], rep(cf$debiased.error[n / 2 + 1:(n / 2)], 10))
  cf.clust$excess.error <- c(cf$excess.error[1:(n / 2)], rep(cf$excess.error[n / 2 + 1:(n / 2)], 10))
  cf.clust$clusters <- c(1:(n / 2), rep(n / 2 + 1:(n / 2), 10))
  cf.clust$equalize.cluster.weights <- TRUE
  tc.clust <- test_calibration(cf.clust)

  expect_equal(tc[, 1:3], tc.clust[, 1:3])
})

test_that("best linear projection is reasonable", {
  n <- 2000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  S <- rbinom(n, 1, 0.5)
  Y <- S * pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  seed = runif(1, 0, .Machine$integer.max)
  forest <- causal_forest(X, Y, W, num.trees = 300, seed = seed)
  forest.shifted.W <- causal_forest(X, Y, W + 1, num.trees = 300, seed = seed)
  forest.2W <- causal_forest(X, Y, W * 2, num.trees = 300, seed = seed)

  # hard-code true regression coefficients
  beta0.true <- 0.2 # approximation to mean(pmax(Z, 0))/2
  beta1.true <- 0.25 # cov(pmax(Z, 0), Z)/2

  blp.all <- best_linear_projection(forest, X[,1:2])
  expect_equal(blp.all[1,1], beta0.true, tolerance = 0.1)
  expect_equal(blp.all[2,1], beta1.true, tolerance = 0.1)
  expect_equal(blp.all[3,1], 0, tolerance = 0.1)

  blp.subset <- best_linear_projection(forest, X[,1:2], subset = (S == 1))
  expect_equal(blp.subset[1,1], beta0.true * 2, tolerance = 0.2)
  expect_equal(blp.subset[2,1], beta1.true * 2, tolerance = 0.2)
  expect_equal(blp.subset[3,1], 0, tolerance = 0.2)

  std.errs <- c((blp.all[1,1] - beta0.true) / blp.all[1,2],
               (blp.all[2,1] - beta1.true) / blp.all[2,2],
               (blp.all[3,1] - 0) / blp.all[3,2],
               (blp.subset[1,1] - beta0.true * 2) / blp.subset[1,2],
               (blp.subset[2,1] - beta1.true * 2) / blp.subset[2,2],
               (blp.subset[3,1] - 0) / blp.subset[3,2])

  expect_lt(mean(abs(std.errs)), 2.5)
  expect_gt(mean(abs(std.errs)), 0.25)

  blp.shifted <- best_linear_projection(forest.shifted.W, X[,1:2])
  expect_equal(blp.all[, "Estimate"], blp.shifted[, "Estimate"], tolerance = 0.05)
  expect_equal(blp.all[, "Std. Error"], blp.shifted[, "Std. Error"], tolerance = 0.05)

  blp.2W <- best_linear_projection(forest.2W, X[,1:2])
  expect_equal(blp.all[, "Estimate"]/2, blp.2W[, "Estimate"], tolerance = 0.05)

  # should be ~ equal to an instrumental forest with Z = W.
  iv.forest <- instrumental_forest(X, Y, W, W,
                                   Y.hat = forest$Y.hat,
                                   W.hat = forest$W.hat,
                                   Z.hat = forest$W.hat,
                                   num.trees = 300,
                                   seed = seed)
  blp.iv.all <- best_linear_projection(iv.forest, X[, 1:2], compliance.score = rep(1, n))
  expect_equal(blp.iv.all[, "Estimate"], blp.all[, "Estimate"], tolerance = 1e-10)
  expect_equal(blp.iv.all[, "Std. Error"], blp.all[, "Std. Error"], tolerance = 1e-10)

  blp.iv.subset <- best_linear_projection(iv.forest, X[, 1:2], compliance.score = rep(1, n), subset = (S == 1))
  expect_equal(blp.iv.subset[, "Estimate"], blp.subset[, "Estimate"], tolerance = 1e-10)
  expect_equal(blp.iv.subset[, "Std. Error"], blp.subset[, "Std. Error"], tolerance = 1e-10)
})

test_that("best linear projection works with edge case input types", {
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  forest <- causal_forest(X, Y, W, num.trees = 50)
  forest.clustered <- causal_forest(X, Y, W, num.trees = 50, clusters = c(rep(1, 100), rep(2, 100)))

  # a single covariate
  blp.single.covariate <- best_linear_projection(forest, X[, 1])
  # a forest trained with clusters, and subset equal to only one of the clusters
  expect_error(best_linear_projection(forest.clustered, X[, 1], subset = which(forest.clustered$clusters == 1)),
               regexp = "The specified subset must contain units from more than one cluster.")
})

test_that("best linear projection works as expected with causal survival forest", {
  n <- 500
  p <- 5
  data <- generate_causal_survival_data(n, p, n.mc = 1, dgp = "simple1")
  data.test <- generate_causal_survival_data(5000, p, n.mc = 10000, dgp = "simple1")
  cs.forest <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max, num.trees = 500)

  ate.true <- mean(data.test$cate)
  blp.ate <- best_linear_projection(cs.forest)
  expect_equal(ate.true, blp.ate[1, "Estimate"], tolerance = 2.1 * sqrt(blp.ate[1, "Std. Error"]))

  A1 <- data.test$X[, 1]
  blp.X1.true <- lm(data.test$cate ~ A1)$coefficients
  blp.X1 <- best_linear_projection(cs.forest, data$X[, 1])
  expect_equal(blp.X1.true[[1]], blp.X1[1, "Estimate"], tolerance = 2.1 * sqrt(blp.X1[1, "Std. Error"]))
  expect_equal(blp.X1.true[[2]], blp.X1[2, "Estimate"], tolerance = 2.1 * sqrt(blp.X1[2, "Std. Error"]))

  weights <- rep(1, n)
  dup <- sample(1:n, 50)
  weights[dup] <- 2
  cs.forest.weight <- causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max,
                                             num.trees = 500, sample.weights = weights)
  cs.forest.dup <- causal_survival_forest(
    rbind(data$X, data$X[dup, ]), c(data$Y, data$Y[dup]), c(data$W, data$W[dup]), c(data$D, data$D[dup]),
    horizon = data$Y.max, num.trees = 500
  )

  blp.weight <- best_linear_projection(cs.forest.weight)
  blp.dup <- best_linear_projection(cs.forest.dup)
  expect_equal(blp.weight[1, "Estimate"], blp.dup[1, "Estimate"], tolerance = 0.01)
})

test_that("best linear projection works as expected with instrumental forest", {
  p <- 5
  n <- 1000
  X <- matrix(2 * runif(n * p) - 1, n, p)
  A <- rnorm(n)
  Z <- rbinom(n, 1, 0.5)
  W <- A + Z * (1 + (X[,2] > 0))
  tau <- X[,1] > 0
  Y <- 2 * (X[,1] <= 0) * A + tau * W + (1 + (sqrt(3) - 1) * (X[,1] > 0)) * rnorm(n)

  iv.forest <- instrumental_forest(X, Y, W, Z, num.trees = 500)
  blp <- best_linear_projection(iv.forest)
  ate <- average_treatment_effect(iv.forest)

  expect_equal(blp[, "Estimate"], ate[["estimate"]], tolerance = 1e-10)
  expect_equal(blp[, "Std. Error"], ate[["std.err"]], tolerance = 1e-4)
})

test_that("ground truth does better than forest-estimates of overlap weights", {
  # A hard setup
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * (p)), n, p)
  eX <- 1 / (1 + exp(-10 * X[, 2]))
  W <- rbinom(n, 1, eX)
  M <- X[, 2]
  TAU <- (1 + X[, 2])^2
  Y <- M + (W - 0.5) * TAU + rnorm(n)

  W.forest <- regression_forest(X, W)
  W.hat.forest <- predict(W.forest)$predictions
  suppressWarnings(W.lm <- glm(W ~ X, family = "binomial"))
  W.hat.lm <- predict(W.lm, type = "response")

  blp.true <- lm(TAU ~ X, weights = eX * (1 - eX))
  blp.lm.wts <- lm(TAU ~ X, weights = W.hat.lm * (1 - W.hat.lm))
  blp.forest.wts <- lm(TAU ~ X, weights = W.hat.forest * (1 - W.hat.forest))

  expect_equal(coef(blp.lm.wts), coef(blp.true), tolerance = 0.1)
  expect_equal(coef(blp.forest.wts), coef(blp.true), tolerance = 0.5)
})

test_that("overlap weighted best linear projection works as expected in a simple setup", {
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * (p)), n, p)
  eX <- 1 / (1 + exp(-3 * X[, 1]))
  W <- rbinom(n, 1, eX)
  TAU <- X[, 2]
  Y <- X[, 1] + (W - 0.5) * TAU + rnorm(n)
  true.blp <- lm(TAU ~ X, weights = eX * (1 - eX))

  cf <- causal_forest(X, Y, W)
  blp.wate <- best_linear_projection(cf, X, target.sample = "overlap")
  expect_equal(blp.wate[3, 1], coef(true.blp)[[3]], tolerance = 3 * blp.wate[3, 2])

  # The overlap-weighted SEs are on average smaller:
  # https://github.com/grf-labs/grf/pull/1258#discussion_r1137992471

  # suppressWarnings(blp <- best_linear_projection(cf, X))
  # se.ratio <- blp.wate[, "Std. Error"] / blp[, "Std. Error"]
  # expect_true(all(se.ratio < 0.95))
})

test_that("causal forest overlap weighted BLP ~same as complete data causal survival BLP", {
  n <- 1000
  p <- 5
  X <- matrix(runif(n * (p)), n, p)
  eX <- 1 / (1 + exp(-3 * X[, 1]))
  W <- rbinom(n, 1, eX)
  TAU <- X[, 2]
  Y <- X[, 1] + W * TAU + runif(n)

  cf <- causal_forest(X, Y, W)
  csf <- causal_survival_forest(X, Y, W, rep(1, n), W.hat = cf$W.hat, horizon = max(Y))

  blp.cf <- best_linear_projection(cf, X, target.sample = "overlap")
  blp.csf <- best_linear_projection(csf, X, target.sample = "overlap")

  expect_equal(blp.cf[, "Estimate"], blp.csf[, "Estimate"], tolerance = 0.05)
  expect_equal(blp.cf[, "Std. Error"], blp.csf[, "Std. Error"], tolerance = 0.05)
})
