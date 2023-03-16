library(grf)

set.seed(1000)

test_that("average effects are translation invariant", {
  p <- 6
  n <- 200
  X <- matrix(2 * runif(n * p) - 1, n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- (X[, 1] > 0) * (2 * W - 1) + 2 * rnorm(n)
  Y.plus.1 <- Y + 1
  forest.causal <- causal_forest(X, Y, W, num.trees = 200, ci.group.size = 4)
  forest.causal.plus.1 <- forest.causal
  forest.causal.plus.1$Y.orig <- forest.causal$Y.orig + 1
  forest.causal.plus.1$Y.hat <- forest.causal$Y.hat + 1

  cate.aipw <- average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  cate.plus.1.aipw <- average_treatment_effect(forest.causal.plus.1, target.sample = "all", method = "AIPW")
  expect_equal(cate.aipw, cate.plus.1.aipw)

  cate.tmle <- average_treatment_effect(forest.causal, target.sample = "all", method = "TMLE")
  cate.plus.1.tmle <- average_treatment_effect(forest.causal.plus.1, target.sample = "all", method = "TMLE")
  expect_equal(cate.tmle, cate.plus.1.tmle)

  catt.aipw <- average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  catt.plus.1.aipw <- average_treatment_effect(forest.causal.plus.1, target.sample = "treated", method = "AIPW")
  expect_equal(catt.aipw, catt.plus.1.aipw)

  catt.tmle <- average_treatment_effect(forest.causal, target.sample = "treated", method = "TMLE")
  catt.plus.1.tmle <- average_treatment_effect(forest.causal.plus.1, target.sample = "treated", method = "TMLE")
  expect_equal(catt.tmle, catt.plus.1.tmle)

  catc.aipw <- average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  catc.plus.1.aipw <- average_treatment_effect(forest.causal.plus.1, target.sample = "control", method = "AIPW")
  expect_equal(catc.aipw, catc.plus.1.aipw)

  catc.tmle <- average_treatment_effect(forest.causal, target.sample = "control", method = "TMLE")
  catc.plus.1.tmle <- average_treatment_effect(forest.causal.plus.1, target.sample = "control", method = "TMLE")
  expect_equal(catc.tmle, catc.plus.1.tmle)

  wate <- average_treatment_effect(forest.causal, target.sample = "overlap")
  wate.plus.1 <- average_treatment_effect(forest.causal.plus.1, target.sample = "overlap")
  expect_equal(wate[1], wate.plus.1[1], tolerance = 0.006)

  # Now test this in "average partial effect" mode. Add some fuzz to the treatments
  # so the "continuous treatment" code path is triggered. However, first cache
  # the debiasing weights so we don't actually have to re-compute them
  debiasing.weights <- (forest.causal$W.orig - forest.causal$W.hat) /
    (forest.causal$W.hat * (1 - forest.causal$W.hat))

  fuzz <- 0.00000001 * rnorm(length(forest.causal$W.orig))
  forest.causal$W.orig <- forest.causal$W.orig + fuzz
  forest.causal.plus.1$W.orig <- forest.causal.plus.1$W.orig + fuzz
  cape <- average_treatment_effect(forest.causal, debiasing.weights=debiasing.weights)
  cape.plus.1 <- average_treatment_effect(forest.causal.plus.1, debiasing.weights=debiasing.weights)
  expect_equal(cape[1], cate.aipw[1])
  expect_equal(cape[1], cape.plus.1[1])
})

test_that("average treatment effect estimates are reasonable", {
  p <- 6
  n <- 2000
  X <- matrix(2 * runif(n * p) - 1, n, p)
  eX <- 0.25 + 0.5 * (X[, 1] > 0)
  W <- rbinom(n, 1, eX)
  TAU <- 4 * (X[, 1] > 0)
  Y <- TAU * (W - 0.5) + rnorm(n)

  forest.causal <- causal_forest(X, Y, W, num.trees = 500, ci.group.size = 1)

  cate.aipw <- average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  expect_equal(cate.aipw[[1]], mean(TAU), tolerance = 0.2)
  expect_equal(cate.aipw[[1]], mean(TAU), tolerance = 3 * cate.aipw[2])

  cate.tmle <- average_treatment_effect(forest.causal, target.sample = "all", method = "TMLE")
  expect_equal(cate.tmle[[1]], mean(TAU), tolerance = 0.2)
  expect_equal(cate.tmle[[1]], mean(TAU), tolerance = 3 * cate.tmle[2])

  expect_equal(cate.aipw[1], cate.tmle[1], tolerance = 0.01)
  expect_equal(cate.aipw[2], cate.tmle[2], tolerance = 0.01)

  catt.aipw <- average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  expect_equal(catt.aipw[[1]], mean(TAU[W == 1]), tolerance = 0.2)
  expect_equal(catt.aipw[[1]], mean(TAU[W == 1]), tolerance = 3 * catt.aipw[2])

  catt.tmle <- average_treatment_effect(forest.causal, target.sample = "treated", method = "TMLE")
  expect_equal(catt.tmle[[1]], mean(TAU[W == 1]), tolerance = 0.2)
  expect_equal(catt.tmle[[1]], mean(TAU[W == 1]), tolerance = 3 * catt.tmle[2])

  expect_equal(catt.aipw[1], catt.tmle[1], tolerance = 0.05)
  expect_equal(catt.aipw[2], catt.tmle[2], tolerance = 0.05)

  catc.aipw <- average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  expect_equal(catc.aipw[[1]], mean(TAU[W == 0]), tolerance = 0.25)
  expect_equal(catc.aipw[[1]], mean(TAU[W == 0]), tolerance = 3 * catc.aipw[2])

  catc.tmle <- average_treatment_effect(forest.causal, target.sample = "control", method = "TMLE")
  expect_equal(catc.tmle[[1]], mean(TAU[W == 0]), tolerance = 0.25)
  expect_equal(catc.tmle[[1]], mean(TAU[W == 0]), tolerance = 3 * catc.tmle[2])

  expect_equal(catc.aipw[1], catc.tmle[1], tolerance = 0.05)
  expect_equal(catc.aipw[2], catc.tmle[2], tolerance = 0.05)

  cape.nocal <- average_treatment_effect(forest.causal)
  expect_equal(cape.nocal[[1]], mean(TAU), tolerance = 0.2)
  expect_equal(cape.nocal[[1]], mean(TAU), tolerance = 3 * cape.nocal[2])
  expect_equal(cate.aipw[1], cape.nocal[1], tolerance = 0.05)
  expect_equal(cate.aipw[2], cape.nocal[2], tolerance = 0.05)

  # The calibration option was eliminated after version 1.2.0.
  # cape.cal <- average_treatment_effect(forest.causal, calibrate.weights = TRUE)
  # expect_true(abs(cape.cal[1] - mean(TAU)) <= 0.2)
  # expect_true(abs(cape.cal[1] - mean(TAU)) <= 3 * cape.cal[2])
  # expect_true(abs(cate.aipw[1] - cape.cal[1]) <= 0.05)
  # expect_true(abs(cate.aipw[2] - cape.cal[2]) <= 0.05)

  wate <- average_treatment_effect(forest.causal, target.sample = "overlap")
  tau.overlap <- sum(eX * (1 - eX) * TAU) / sum(eX * (1 - eX))
  expect_equal(wate[[1]], tau.overlap, tolerance = 0.2)
  expect_equal(wate[[1]], tau.overlap, tolerance = 3 * wate[2])

  cate.aipw.pos <- average_treatment_effect(forest.causal,
    target.sample = "all",
    method = "AIPW", subset = X[, 1] > 0
  )
  cate.tmle.pos <- average_treatment_effect(forest.causal,
    target.sample = "all",
    method = "TMLE", subset = X[, 1] > 0
  )
  cate.aipw.pos.treat <- average_treatment_effect(forest.causal,
    target.sample = "treated",
    method = "AIPW", subset = which(X[, 1] > 0)
  )
  cate.tmle.pos.control <- average_treatment_effect(forest.causal,
    target.sample = "control",
    method = "TMLE", subset = which(X[, 1] > 0)
  )
  wate.pos <- average_treatment_effect(forest.causal, target.sample = "overlap", subset = X[, 1] > 0)

  expect_equal(cate.aipw.pos[[1]], 4, tolerance = 0.2)
  expect_equal(cate.tmle.pos[[1]], 4, tolerance = 0.2)
  expect_equal(cate.aipw.pos.treat[[1]], 4, tolerance = 0.2)
  expect_equal(cate.tmle.pos.control[[1]], 4, tolerance = 0.3)
  expect_equal(wate.pos[[1]], 4, tolerance = 0.2)
})

test_that("average partial effect estimates are reasonable", {
  p <- 6
  n <- 1000
  X <- matrix(2 * runif(n * p) - 1, n, p)
  eX <- 0.25 + 0.5 * (X[, 1] > 0)
  W <- rbinom(n, 1, eX) + rnorm(n)
  TAU <- 4 * (X[, 1] > 0)
  Y <- TAU * (W - 0.5) + rnorm(n)
  forest.causal <- causal_forest(X, Y, W,
    num.trees = 500,
    ci.group.size = 1, clusters = rep(1:(n / 2), 2)
  )
  cape.pos <- average_treatment_effect(forest.causal, subset = X[, 1] > 0)
  expect_equal(cape.pos[["estimate"]], 4, tolerance = 0.1)
})

test_that("average treatment effects larger example works", {
  n <- 4000
  p <- 10

  X <- matrix(runif(n * p), n, p)
  E <- (0.4 + dbeta(X[, 2], 2, 4)) / 4
  W <- rbinom(n, 1, E)
  M <- 2 * X[, 2] - 1
  TAU <- (1 + 1 / (1 + exp(-20 * (X[, 1] - 0.3)))) * (1 + 1 / (1 + exp(-20 * (X[, 2] - 0.3))))
  Y <- M + (W - 0.5) * TAU + rnorm(n)

  forest.causal <- causal_forest(X, Y, W, num.trees = 1000, ci.group.size = 1)

  cate.aipw <- average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  expect_equal(cate.aipw[[1]], mean(TAU), tolerance = 3 * cate.aipw[2])

  cate.tmle <- average_treatment_effect(forest.causal, target.sample = "all", method = "TMLE")
  expect_equal(cate.tmle[[1]], mean(TAU), tolerance = 3 * cate.tmle[2])

  catt.aipw <- average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  expect_equal(catt.aipw[[1]], mean(TAU[W == 1]), tolerance = 3 * catt.aipw[2])

  catt.tmle <- average_treatment_effect(forest.causal, target.sample = "treated", method = "TMLE")
  expect_equal(catt.tmle[[1]], mean(TAU[W == 1]), tolerance = 3 * catt.tmle[2])

  catc.aipw <- average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  expect_equal(catc.aipw[[1]], mean(TAU[W == 0]), tolerance = 3 * catc.aipw[2])

  catc.tmle <- average_treatment_effect(forest.causal, target.sample = "control", method = "TMLE")
  expect_equal(catc.tmle[[1]], mean(TAU[W == 0]), tolerance = 3 * catc.tmle[2])
})

test_that("average partial effects larger example works", {
  n <- 4000
  p <- 10

  X <- matrix(runif(n * p), n, p)
  E <- (0.4 + dbeta(X[, 2], 2, 4)) / 4
  W <- rbinom(n, 1, E) + 0.2 * rnorm(n)
  M <- 2 * X[, 2] - 1
  TAU <- (1 + 1 / (1 + exp(-20 * (X[, 1] - 0.3)))) * (1 + 1 / (1 + exp(-20 * (X[, 2] - 0.3))))
  Y <- M + (W - 0.5) * TAU + rnorm(n)

  forest.causal <- causal_forest(X, Y, W, num.trees = 1000, ci.group.size = 1)

  cape <- average_treatment_effect(forest.causal)
  expect_equal(cape[[1]],  mean(TAU), tolerance = 0.2)
  expect_equal(cape[[1]], mean(TAU), tolerance = 3 * cape[2])
})

test_that("average treatment effect with overlap: larger example works", {
  n <- 4000
  p <- 10

  X <- matrix(rnorm(n * (p)), n, p)
  eX <- 1 / (1 + exp(-10 * X[, 2]))
  W <- rbinom(n, 1, eX)
  M <- X[, 2]
  TAU <- (1 + X[, 2])^2
  Y <- M + (W - 0.5) * TAU + rnorm(n)

  forest.causal <- causal_forest(X, Y, W, num.trees = 1000, ci.group.size = 1)

  wate <- average_treatment_effect(forest.causal, target.sample = "overlap")
  tau.overlap <- sum(eX * (1 - eX) * TAU) / sum(eX * (1 - eX))
  expect_equal(wate[[1]], tau.overlap, tolerance = 0.2)
  expect_equal(wate[[1]], tau.overlap, tolerance = 3 * wate[2])

  # Should be very similar to a best linear projection on a constant with target.sample = "overlap"
  blp.wate <- best_linear_projection(forest.causal, target.sample = "overlap")
  expect_equal(blp.wate[1, 1], wate[[1]], tolerance = 0.05)
  expect_equal(blp.wate[1, 2], wate[[2]], tolerance = 0.05)
})

test_that("cluster robust average effects are consistent", {

  p <- 3
  n <- 800

  X <- matrix(2 * runif(n * p) - 1, n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- (X[, 1] > 0) * (2 * W - 1) + 2 * rnorm(n)

  Xc <- rbind(X, X, X, X)
  Wc <- c(W, W, W, W)
  Yc <- c(Y, Y, Y, Y)
  clust <- c(1:n, 1:n, 1:n, 1:n)

  forest.causal <- causal_forest(X, Y, W, num.trees = 2000, ci.group.size = 4)
  forest.causal.clust <- causal_forest(Xc, Yc, Wc,
                                      num.trees = 2000,
                                      ci.group.size = 4, clusters = clust)

  cate.aipw <- average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  cate.clust.aipw <- average_treatment_effect(forest.causal.clust, target.sample = "all", method = "AIPW")
  expect_equal(cate.aipw[1], cate.clust.aipw[1], tolerance = 0.05)
  expect_equal(cate.aipw[2], cate.clust.aipw[2], tolerance = 0.008)

  catt.aipw <- average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  catt.clust.aipw <- average_treatment_effect(forest.causal.clust, target.sample = "treated", method = "AIPW")
  expect_equal(catt.aipw[1], catt.clust.aipw[1], tolerance = 0.05)
  expect_equal(catt.aipw[2], catt.clust.aipw[2], tolerance = 0.008)

  catc.aipw <- average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  catc.clust.aipw <- average_treatment_effect(forest.causal.clust, target.sample = "control", method = "AIPW")
  expect_equal(catc.aipw[1], catc.clust.aipw[1], tolerance = 0.05)
  expect_equal(catc.aipw[2], catc.clust.aipw[2], tolerance = 0.008)

  cape <- average_treatment_effect(forest.causal, num.trees.for.weights = 200)
  cape.clust <- average_treatment_effect(forest.causal.clust, num.trees.for.weights = 200)
  expect_equal(cape[1], cape.clust[1], tolerance = 0.05)
  expect_equal(cape[2], cape.clust[2], tolerance = 0.005)

  wate <- average_treatment_effect(forest.causal, target.sample = "overlap")
  wate.clust <- average_treatment_effect(forest.causal.clust, target.sample = "overlap")
  expect_equal(wate[1], wate.clust[1], tolerance = 0.05)
  expect_equal(wate[2], wate.clust[2], tolerance = 0.005)
})

test_that("cluster robust average effects do weighting correctly", {
  t0 <- 2
  K <- 400
  n <- 11 * K
  p <- 4
  clust <- 1:n %% K + K * as.numeric(1:n >= K)
  tau <- 2 * t0 * as.numeric(clust < K)

  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- tau * W + 2 * rnorm(n)

  forest.causal <- causal_forest(X, Y, W, clusters = clust, equalize.cluster.weights = TRUE, num.trees = 400)

  cate.aipw <- average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  expect_lte(abs(cate.aipw[1] - t0) / (3 * cate.aipw[2]),  1)
  expect_lte(cate.aipw[2], 0.2)

  # The best linear projection with NULL covariates should match the ATE estimate via AIPW.
  # The reason the standard error estimates don't match exactly is that the function
  # `best_linear_projection` estimates standard errors using the more general function
  # `coeftest`, whereas `average_treatment_effect` uses a direct calculation.
  cate.aipw.blp <- best_linear_projection(forest.causal, A = NULL)
  expect_equal(as.numeric(cate.aipw[1]), cate.aipw.blp[1,1])
  expect_equal(as.numeric(cate.aipw[2]), cate.aipw.blp[1,2], tolerance = 0.0001)

  catt.aipw <- average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  expect_lte(abs(catt.aipw[1] - t0) / (3 * catt.aipw[2]), 1)
  expect_lte(catt.aipw[2], 0.2)

  catc.aipw <- average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  expect_lte(abs(catc.aipw[1] - t0) / (3 * catc.aipw[2]), 1)
  expect_lte(catc.aipw[2], 0.2)

  cape <- average_treatment_effect(forest.causal, num.trees.for.weights = 200)
  expect_lte(abs(cape[1] - t0) / (3 * cape[2]), 1)
  expect_lte(cape[2], 0.2)

  wate <- average_treatment_effect(forest.causal, target.sample = "overlap")
  expect_lte(abs(wate[1] - t0) / (3 * wate[2]), 1)
  expect_lte(wate[2], 0.2)

  # An IV forest with W = Z should behave just like a causal forest
  # with treatment W.
  forest.instrumental <- instrumental_forest(X, Y, W, W,
                                             Y.hat = forest.causal$Y.hat,
                                             W.hat = forest.causal$W.hat,
                                             Z.hat = forest.causal$W.hat,
                                             clusters = clust,
                                             equalize.cluster.weights = TRUE,
                                             num.trees = 400)
  compliance.score <- rep(1, n)
  aclate <- average_treatment_effect(forest.instrumental, compliance.score=compliance.score)
  expect_equal(cate.aipw["estimate"], aclate["estimate"], tolerance = 0.03)
  expect_equal(cate.aipw["std.err"], aclate["std.err"], tolerance = 0.002)
})

test_that("cluster robust average effects do weighting correctly with IPCC weights", {
  # This test may warn about low overlap, so we disable warnings here to avoid
  # test failures, and restore it at the end.
  options(warn = -1)

  t0 <- 2
  t1 <- 3
  K <- 100
  p <- 4
  cluster.sizes <- pmax(20, round(40 + 40 * rt(K, df = 3)))
  n <- sum(cluster.sizes)
  clust <- rep(1:K, cluster.sizes)

  X <- matrix(rnorm(n * p), n, p)
  e <- 1 / (1 + exp(-X[, 1] + X[, 2]))
  W <- rbinom(n, 1, e)
  tau <- 2 * t0 * as.numeric(cluster.sizes[clust] >= median(cluster.sizes)) +
    2 * t1 * as.numeric(X[, 3] > 0)
  Y <- 2 * X[, 1] + tau * W + 2 * rnorm(n)

  # strictly speaking, we should add variance when comparing to this
  # because the variances we report are for the conditional average and this is the population average.
  true.ate <- t0 + t1

  e.cc <- 0.2 + 0.6 * as.numeric(X[, 3] > 0)
  cc <- as.logical(rbinom(n, 1, e.cc))
  ipcc.sample.weights <- 1 / e.cc
  clust.sample.weights <- 1 / cluster.sizes[clust]
  sample.weights <- ipcc.sample.weights * clust.sample.weights

  forest.weighted <- causal_forest(
    X[cc, ], Y[cc], W[cc],
    sample.weights = sample.weights[cc],
    clusters = clust[cc], equalize.cluster.weights = FALSE, num.trees = 400
  )
  forest.unweighted <- causal_forest(X[cc, ], Y[cc], W[cc], clusters = clust[cc],
                                     equalize.cluster.weights = FALSE,  num.trees = 400)

  cate.aipw <- average_treatment_effect(forest.weighted, target.sample = "all", method = "AIPW")
  biased.cate.aipw <- average_treatment_effect(forest.unweighted, target.sample = "all", method = "AIPW")
  expect_lte(abs(cate.aipw[1] - true.ate) / (3 * cate.aipw[2]), 1)
  expect_gte(abs(biased.cate.aipw[1] - true.ate) / (3 * biased.cate.aipw[2]), 1)

  # The best linear projection with NULL covariates should match the ATE estimate via AIPW.
  # The reason the standard error estimates don't match exactly is that the function
  # `best_linear_projection` estimates standard errors using the more general function
  # `coeftest`, whereas `average_treatment_effect` uses a direct calculation.
  cate.aipw.blp <- best_linear_projection(forest.weighted, A = NULL)
  expect_equal(as.numeric(cate.aipw[1]), cate.aipw.blp[1,1])
  expect_equal(as.numeric(cate.aipw[2]), cate.aipw.blp[1,2], tolerance = 0.01)
  biased.cate.aipw.blp <- best_linear_projection(forest.unweighted, A = NULL)
  expect_equal(as.numeric(biased.cate.aipw[1]), biased.cate.aipw.blp[1,1])
  expect_equal(as.numeric(biased.cate.aipw[2]), biased.cate.aipw.blp[1,2], tolerance = 0.01)

  catt.aipw <- average_treatment_effect(forest.weighted, target.sample = "treated", method = "AIPW")
  biased.catt.aipw <- average_treatment_effect(forest.unweighted, target.sample = "treated", method = "AIPW")
  expect_lte(abs(catt.aipw[1] - true.ate) / (3 * catt.aipw[2]), 1)
  expect_gte(abs(biased.cate.aipw[1] - true.ate) / (3 * biased.catt.aipw[2]), 1)

  catc.aipw <- average_treatment_effect(forest.weighted, target.sample = "control", method = "AIPW")
  biased.catc.aipw <- average_treatment_effect(forest.unweighted, target.sample = "control", method = "AIPW")
  expect_lte(abs(catc.aipw[1] - true.ate) / (3 * catc.aipw[2]), 1)
  expect_gte(abs(biased.catc.aipw[1] - true.ate) / (3 * biased.catc.aipw[2]), 1)

  cape <- average_treatment_effect(forest.weighted, num.trees.for.weights = 200)
  biased.cape <- average_treatment_effect(forest.unweighted, num.trees.for.weights = 200)
  expect_lte(abs(cape[1] - true.ate) / (3 * cape[2]), 1)
  expect_gte(abs(biased.cape[1] - true.ate) / (3 * biased.cape[2]), 1)

  wate <- average_treatment_effect(forest.weighted, target.sample = "overlap")
  biased.wate <- average_treatment_effect(forest.unweighted, target.sample = "overlap")
  expect_lte(abs(wate[1] - true.ate) / (3 * wate[2]), 1)
  expect_gte(abs(biased.wate[1] - true.ate) / (3 * biased.wate[2]), 1)

  options(warn = 2)
})

test_that("average effect estimation doesn't error on data with a single feature", {
  p <- 1
  n <- 100

  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)
  W <- rbinom(n, size = 1, prob = 0.5)

  forest <- causal_forest(X, Y, W)
  average_treatment_effect(forest)
  expect_true(TRUE) # so we don't get a warning about an empty test
})

test_that("average conditional local average treatment effect estimation is reasonable", {
  p <- 10
  n <- 1000

  X <- matrix(2 * runif(n * p) - 1, n, p)
  A <- rnorm(n)
  Z <- rbinom(n, 1, 0.5)
  W <- A + Z * (1 + (X[,2] > 0))
  tau <- X[,1] > 0
  Y <- 2 * (X[,1] <= 0) * A + tau * W + (1 + (sqrt(3) - 1) * (X[,1] > 0)) * rnorm(n)

  forest.iv <- instrumental_forest(X, Y, W, Z, num.trees = 250)
  compliance.forest <- causal_forest(forest.iv$X.orig,
                                     Y=forest.iv$W.orig,
                                     W=forest.iv$Z.orig,
                                     Y.hat=forest.iv$W.hat,
                                     W.hat=forest.iv$Z.hat,
                                     sample.weights = forest.iv$sample.weights,
                                     num.trees = 250)
  compliance.score <- predict(compliance.forest)$predictions

  tau.hat <- average_treatment_effect(forest.iv, compliance.score=compliance.score)
  tau.x1p <- average_treatment_effect(forest.iv, compliance.score=compliance.score,
                          subset = X[,1] > 0)

  expect_equal(as.numeric(tau.hat["estimate"]), mean(tau), tolerance = 0.2)
  expect_lt(abs(tau.hat["estimate"] - mean(tau)) / tau.hat["std.err"], 3)

  expect_equal(as.numeric(tau.x1p["estimate"]), mean(tau[X[,1] > 0]), tolerance = 0.3)
  expect_lt(abs(tau.x1p["estimate"] - mean(tau[X[,1] > 0])) / tau.x1p["std.err"], 3)
})

test_that("average effect estimation handles SE's with sample weights=0 consistently", {
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  wts <- sample(c(0, 1, 2), n, replace = TRUE)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, sample.weights = wts)

  expect_equal(average_treatment_effect(cf),
               average_treatment_effect(cf, subset = wts > 0))
  expect_equal(average_treatment_effect(cf, target.sample = "control"),
               average_treatment_effect(cf, target.sample = "control", subset = wts > 0))
  expect_equal(average_treatment_effect(cf, target.sample = "treated"),
               average_treatment_effect(cf, target.sample = "treated", subset = wts > 0))

  cl <- sample(c(5:20), n, replace = TRUE)
  cf.clust <- causal_forest(X, Y, W, W.hat = 0.5, sample.weights = wts, clusters = cl)

  expect_equal(average_treatment_effect(cf.clust),
               average_treatment_effect(cf.clust, subset = wts > 0))
  expect_equal(average_treatment_effect(cf.clust, target.sample = "control"),
               average_treatment_effect(cf.clust, target.sample = "control", subset = wts > 0))
  expect_equal(average_treatment_effect(cf.clust, target.sample = "treated"),
               average_treatment_effect(cf.clust, target.sample = "treated", subset = wts > 0))
})
