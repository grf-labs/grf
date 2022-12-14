test_that("rank_average_treatment_effect works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + pmax(X[, 3], 0) + runif(n)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)

  prio <- get_scores(cf)
  rate <- rank_average_treatment_effect(cf, prio)
  capture.output(print(rate))
  plot(rate)
  plot(rate, ylab = "42")
  plot(rate, ylab = "42", col = "green", sub = "sub")
  plot(rate, xlab = "treated frac.", legend.args = list(legend = "a prio"),
       abline.args = list(lty = 1), ci.args = list(lty = 3))

  q.length <- nrow(rate$TOC)
  expect_equal(rate$TOC[q.length, "estimate"], 0, tolerance = 1e-10) # Last TOC curve entry (q=1) = zero.
  expect_equal(rate$TOC[q.length, "std.err"], 0, tolerance = 1e-10) # Last TOC curve entry = zero.
  expect_true(all(rate$TOC$estimate[1:(q.length - 1)] >= 0)) # for this prio all points on curve are > 0

  rate.lo <- rank_average_treatment_effect(cf, prio, subset = prio < 0)
  expect_gt(rate[["estimate"]], rate.lo[["estimate"]])

  rate.all.eq <- rank_average_treatment_effect(cf, rep(1, n))
  expect_equal(rate.all.eq[["estimate"]], 0, tolerance = 1e-10)

  rand.prio <- sample(1:100, n, TRUE)
  autoc.rand <- rank_average_treatment_effect(cf, rand.prio)
  expect_equal(autoc.rand[["estimate"]], 0, tolerance = 3 * autoc.rand[["std.err"]])

  qini.rand <- rank_average_treatment_effect(cf, rand.prio, target = "QINI")
  expect_equal(qini.rand[["estimate"]], 0, tolerance = 3 * qini.rand[["std.err"]])

  cf.survival <- causal_survival_forest(X, Y, W, rep(1, n), W.hat = 0.5, num.trees = 250, horizon = max(Y))
  autoc.cfs <- rank_average_treatment_effect(cf.survival, rand.prio)
  expect_equal(autoc.cfs[["estimate"]], 0, tolerance = 3 * autoc.cfs[["std.err"]])
})

test_that("rank_average_treatment_effect discriminates between HTE/No HTE", {
  # HTEs
  n <- 1000
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + pmax(X[, 3], 0) + runif(n)
  train <- 1:500
  eval <- -train

  cf.train <- causal_forest(X[train, ], Y[train], W[train], W.hat = 0.5, num.trees = 250)
  prio.cate <- predict(cf.train, X[eval, ])$predictions
  prio.rand <- runif(length(prio.cate))
  cf.eval <- causal_forest(X[eval, ], Y[eval], W[eval], W.hat = 0.5, num.trees = 250)

  rate.hte <- rank_average_treatment_effect(cf.eval, cbind(prio.cate, prio.rand))
  expect_gt(rate.hte$estimate[[1]] / rate.hte$std.err[[1]], 2.5)
  expect_equal(rate.hte$estimate[[2]], 0, tolerance = 2.5 * rate.hte$std.err[[2]])
  expect_gt(rate.hte$estimate[[3]] / rate.hte$std.err[[3]], 2.5)

  # No HTEs
  tau <- 0
  Y <- tau * W + pmax(X[, 3], 0) + runif(n)

  cf.train <- causal_forest(X[train, ], Y[train], W[train], W.hat = 0.5, num.trees = 250)
  prio.cate <- predict(cf.train, X[eval, ])$predictions
  prio.rand <- runif(length(prio.cate))
  cf.eval <- causal_forest(X[eval, ], Y[eval], W[eval], W.hat = 0.5, num.trees = 250)

  rate.nohte <- rank_average_treatment_effect(cf.eval, cbind(prio.cate, prio.rand))
  expect_equal(rate.nohte$estimate[[1]], 0, tolerance = 2.5 * rate.nohte$std.err[[1]])
  expect_equal(rate.nohte$estimate[[2]], 0, tolerance = 2.5 * rate.nohte$std.err[[2]])
  expect_equal(rate.nohte$estimate[[3]], 0, tolerance = 2.5 * rate.nohte$std.err[[3]])
})

test_that("TOC grid works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)
  prio <- tau

  # Computing TOC on grid 1/n <= q <= 1 agrees exactly with AUTOC.
  q.full <- seq(1/n, 1, by = 1/n)
  rate.full <-rank_average_treatment_effect(cf, prio, q = q.full, R = 0)
  autoc <- rate.full$estimate
  expect_equal(mean(rate.full$TOC$estimate), autoc, tolerance = 1e-10)

  # Can ask for reasonable coverage of the entire curve.
  rand.prio <- sample(1:100, n, TRUE)
  rate <- rank_average_treatment_effect(cf, rand.prio)
  TOC <- rate$TOC$estimate
  TOC.se <- rate$TOC$std.err
  for (i in seq_along(TOC)) {
    expect_equal(TOC[i], 0, tolerance = 3.1 * TOC.se[i])
  }

  q5 <- seq(0.05, 1, by = 0.05)
  rate.q5 <- rank_average_treatment_effect(cf, rand.prio, q = q5)
  TOC.q5 <- rate.q5$TOC$estimate
  TOC.q5.se <- rate.q5$TOC$std.err
  for (i in seq_along(TOC.q5)) {
    expect_equal(TOC.q5[i], 0, tolerance = 3.1 * TOC.q5.se[i])
  }
})

test_that("sample weighted TOC grid works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  wts <- sample(c(1, 2), n, TRUE)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250, sample.weights = wts)
  prio <- tau
  sort.order <- order(prio, decreasing = TRUE)

  # Computing TOC on grid 1/n <= q <= 1 agrees exactly with AUTOC.
  q.full <- cumsum(wts[sort.order] / sum(wts))
  q.full[length(q.full)] <- 1 # In case ~eps difference from 1.
  rate.full <-rank_average_treatment_effect(cf, prio, q = q.full, R = 0)
  autoc <- rate.full$estimate
  expect_equal(weighted.mean(rate.full$TOC$estimate, wts[sort.order]), autoc, tolerance = 1e-10)

  # Can ask for reasonable coverage of the entire curve.
  rand.prio <- sample(1:100, n, TRUE)
  rate <- rank_average_treatment_effect(cf, rand.prio)
  TOC <- rate$TOC$estimate
  TOC.se <- rate$TOC$std.err
  for (i in seq_along(TOC)) {
    expect_equal(TOC[i], 0, tolerance = 3.1 * TOC.se[i])
  }

  q5 <- seq(0.05, 1, by = 0.05)
  rate.q5 <- rank_average_treatment_effect(cf, rand.prio, q = q5)
  TOC.q5 <- rate.q5$TOC$estimate
  TOC.q5.se <- rate.q5$TOC$std.err
  for (i in seq_along(TOC.q5)) {
    expect_equal(TOC.q5[i], 0, tolerance = 3.1 * TOC.q5.se[i])
  }
})

test_that("rank_average_treatment_effect agrees with plain brute-force calculation", {
  n <- 50
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)
  DR.scores <- get_scores(cf)

  # 1. Unique priorities
  prio <- runif(n)
  autoc <- rank_average_treatment_effect(cf, prio, R = 0)
  qini <- rank_average_treatment_effect(cf, prio, target = "QINI", R = 0)

  sort.idx <- order(prio, decreasing = TRUE)
  TOC <- rep(NA, n)
  for (i in 1:n) {
    TOC[i] <- mean(DR.scores[sort.idx[1:i]]) - mean(DR.scores)
  }
  AUTOC <- mean(TOC)
  QINI <- sum(seq.int(1, length(TOC)) / length(TOC)^2 * TOC)

  expect_equal(autoc[["estimate"]], AUTOC, tolerance = 1e-10)
  expect_equal(qini[["estimate"]], QINI, tolerance = 1e-10)

  # 2. Duplicate priorities
  prio.dup <- sample(c(13, 14, 15, 18, 21, 58), n, replace = TRUE)
  autoc.dup <- rank_average_treatment_effect(cf, prio.dup, R = 0)
  qini.dup <- rank_average_treatment_effect(cf, prio.dup, target = "QINI", R = 0)

  # average the doubly robust scores within tied groups
  scores.by.prio <- split(DR.scores, prio.dup) # orders by prio in increasing order
  ties.count.by.prio <- lapply(scores.by.prio, length)
  mean.scores.by.prio <- lapply(scores.by.prio, mean)
  # scores in decreasing priority order
  scores.order <- rev(unlist(rep(mean.scores.by.prio, ties.count.by.prio)))

  TOC.dup <- rep(NA, n)
  for (i in 1:n) {
    TOC.dup[i] <- mean(scores.order[1:i]) - mean(DR.scores)
  }
  AUTOC.dup <- mean(TOC.dup)
  QINI.dup <- sum(seq.int(1, length(TOC.dup)) / length(TOC.dup)^2 * TOC.dup)

  expect_equal(autoc.dup[["estimate"]], AUTOC.dup, tolerance = 1e-10)
  expect_equal(qini.dup[["estimate"]], QINI.dup, tolerance = 1e-10)
})

test_that("sample weighted rank_average_treatment_effect agrees with plain brute-force calculation", {
  n <- 50
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  wts <- sample(c(1, 2), n, TRUE)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250, sample.weights = wts)
  DR.scores <- get_scores(cf)

  # 1. Unique priorities
  prio <- runif(n)
  autoc <- rank_average_treatment_effect(cf, prio, R = 0)
  qini <- rank_average_treatment_effect(cf, prio, target = "QINI", R = 0)

  sort.idx <- order(prio, decreasing = TRUE)
  TOC <- rep(NA, n)
  for (i in 1:n) {
    TOC[i] <- weighted.mean(DR.scores[sort.idx[1:i]], wts[sort.idx[1:i]]) - weighted.mean(DR.scores, wts)
  }
  AUTOC <- weighted.mean(TOC, wts[sort.idx])
  QINI <- sum(cumsum(wts[sort.idx]) / sum(wts) * wts[sort.idx] * TOC) / sum(wts)

  expect_equal(autoc[["estimate"]], AUTOC, tolerance = 1e-10)
  expect_equal(qini[["estimate"]], QINI, tolerance = 1e-10)

  # 2. Duplicate priorities
  prio.dup <- sample(c(13, 14, 15, 18, 21, 58), n, replace = TRUE)
  autoc.dup <- rank_average_treatment_effect(cf, prio.dup, R = 0)
  qini.dup <- rank_average_treatment_effect(cf, prio.dup, target = "QINI", R = 0)

  # sample weights-average the doubly robust scores within tied groups
  sample.by.prio <- split(seq_along(DR.scores), prio.dup) # orders by prio in increasing order
  ties.count.by.prio <- lapply(sample.by.prio, length)
  mean.scores.by.prio <- lapply(sample.by.prio, function(i) weighted.mean(DR.scores[i], wts[i]))
  # scores in decreasing priority order
  scores.order <- rev(unlist(rep(mean.scores.by.prio, ties.count.by.prio)))
  sort.idx.dup <- unlist(rev(sample.by.prio))

  TOC.dup <- rep(NA, n)
  for (i in 1:n) {
    TOC.dup[i] <- weighted.mean(scores.order[1:i], wts[sort.idx.dup[1:i]]) - weighted.mean(DR.scores, wts)
  }
  AUTOC.dup <- weighted.mean(TOC.dup, wts[sort.idx.dup])
  QINI.dup <- sum(cumsum(wts[sort.idx.dup]) / sum(wts) * wts[sort.idx.dup] * TOC.dup) / sum(wts)

  expect_equal(autoc.dup[["estimate"]], AUTOC.dup, tolerance = 1e-10)
  expect_equal(qini.dup[["estimate"]], QINI.dup, tolerance = 1e-10)
})

test_that("sample weighted rank_average_treatment_effect is invariant to weight scaling", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  wts <- sample(c(1, 2), n, TRUE)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250, sample.weights = wts)

  set.seed(42)
  autoc <- rank_average_treatment_effect(cf, tau, target = "AUTOC", R = 10)
  set.seed(42)
  qini <- rank_average_treatment_effect(cf, tau, target = "QINI", R = 10)
  cf$sample.weights = cf$sample.weights / 1000
  set.seed(42)
  autoc.s <- rank_average_treatment_effect(cf, tau, target = "AUTOC", R = 10)
  set.seed(42)
  qini.s <- rank_average_treatment_effect(cf, tau, target = "QINI", R = 10)

  expect_equal(autoc, autoc.s)
  expect_equal(qini, qini.s)
})

test_that("sample weighted rank_average_treatment_effect is approx. same as duplication", {
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

  # Duplicating some samples and giving same samples weight 2 gives very similar AUTOC
  dupe <- sample(1:250, 150)
  X.dupe <- rbind(X, X[dupe,])
  W.dupe <- c(W, W[dupe])
  tau.dupe <- c(tau, tau[dupe])
  Y.dupe <- c(Y, Y[dupe])
  wts <- rep(1, n)
  wts[dupe] <- 2

  cf <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0.5, num.trees = 250, sample.weights = wts)
  cf.dupe <- causal_forest(X.dupe, Y.dupe, W.dupe, Y.hat = 0, W.hat = 0.5, num.trees = 250)

  autoc <- rank_average_treatment_effect(cf, tau, target = "AUTOC", R = 0)
  qini <- rank_average_treatment_effect(cf, tau, target = "QINI", R = 0)
  autoc.dupe <- rank_average_treatment_effect(cf.dupe, tau.dupe, target = "AUTOC", R = 0)
  qini.dupe <- rank_average_treatment_effect(cf.dupe, tau.dupe, target = "QINI", R = 0)

  expect_equal(autoc$estimate, autoc.dupe$estimate, tolerance = 0.01)
  expect_equal(qini$estimate, qini.dupe$estimate, tolerance = 0.01)
  expect_equal(autoc$TOC$estimate, autoc.dupe$TOC$estimate, tolerance = 0.01)
  expect_equal(autoc$TOC[nrow(autoc$TOC), "estimate"], 0, tolerance = 1e-10) # Last TOC curve entry (q=1) = zero.
})

test_that("IPCC weighting rank_average_treatment_effect improves complete-data RATE estimate.", {
  n <- 3000
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  e <- 0.5
  W <- rbinom(n, 1, e)
  tau <- 2 * (X[, 1] > 0 & X[, 5] > 0) -
    0.5 * (X[, 2] > 0) - 0.5 * (X[, 3] > 0) - 0.5 * (X[, 4] > 0)
  Y <- W * tau + rnorm(n)

  e.cc <- 1 - 0.9 * (X[, 1] > 0 & X[, 5] > 0)
  cc <- as.logical(rbinom(n, 1, e.cc))
  sample.weights <- 1 / e.cc

  cf.full <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0.5, num.trees = 250)
  cf.ipc <- causal_forest(X[cc, ], Y[cc], W[cc], Y.hat = 0, W.hat = 0.5, sample.weights = sample.weights[cc], num.trees = 250)
  cf.nowt <- causal_forest(X[cc, ], Y[cc], W[cc], Y.hat = 0, W.hat = 0.5, num.trees = 250)

  rate.full <- rank_average_treatment_effect(cf.full, tau, R = 0)
  rate.ipc <- rank_average_treatment_effect(cf.ipc, tau[cc], R = 0)
  rate.nowt <- rank_average_treatment_effect(cf.nowt, tau[cc], R = 0)

  expect_gt(rate.full$estimate / rate.nowt$estimate, 1.2 * rate.full$estimate / rate.ipc$estimate)
})

test_that("rank_average_treatment_effect is invariant to permuting samples", {
  n <- 10
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  S <- c(10, 10, 10, 9, 8, 8, 7, 7, 7, 6)
  q <- c(0.125, 0.3, 0.375, 0.5, 0.75, 0.8, 1)
  cf <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0.5, num.trees = 250)

  j <- c(2, 3, 1, 4, 6, 5, 8, 9, 7, 10)
  cf.perm <- cf
  cf.perm$Y.orig <- cf$Y.orig[j]
  cf.perm$W.orig <- cf$W.orig[j]
  cf.perm$predictions <- cf$predictions[j, , drop = FALSE]
  cf.perm$sample.weights <- cf$sample.weights[j]
  expect_equal(get_scores(cf.perm), get_scores(cf)[j])

  r1 <- rank_average_treatment_effect(cf, S, q = q)
  r2 <- rank_average_treatment_effect(cf.perm, S, q = q)

  expect_equal(r1$TOC$estimate, r2$TOC$estimate, tolerance = 1e-15)
  expect_equal(r1$estimate, r2$estimate, tolerance = 1e-15)
  expect_equal(r1$TOC$estimate[length(q)], 0, tolerance = 1e-15)
})

test_that("sample weighted rank_average_treatment_effect is invariant to permuting samples", {
  n <- 10
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  S <- c(10, 10, 10, 9, 8, 8, 7, 7, 7, 6)
  wts <- c(1, 3, 1, 1, 0.5, 1.5, 1, 1, 2, 4)
  q <- c(0.125, 0.3, 0.375, 0.5, 0.75, 0.8, 1)
  cf <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0.5, sample.weights = wts, num.trees = 250)

  j <- c(2, 3, 1, 4, 6, 5, 8, 9, 7, 10)
  cf.perm <- cf
  cf.perm$Y.orig <- cf$Y.orig[j]
  cf.perm$W.orig <- cf$W.orig[j]
  cf.perm$predictions <- cf$predictions[j, , drop = FALSE]
  cf.perm$sample.weights <- cf$sample.weights[j]
  expect_equal(get_scores(cf.perm), get_scores(cf)[j])

  r1 <- rank_average_treatment_effect(cf, S, q = q)
  r1.qini <- rank_average_treatment_effect(cf, S, q = q, target = "QINI")
  r2 <- rank_average_treatment_effect(cf.perm, S, q = q)
  r2.qini <- rank_average_treatment_effect(cf.perm, S, q = q, target = "QINI")

  expect_equal(r1$TOC$estimate, r2$TOC$estimate, tolerance = 1e-15)
  expect_equal(r1$estimate, r2$estimate, tolerance = 1e-2) # ascribed to finite precision...
  expect_equal(r1.qini$estimate, r2.qini$estimate, tolerance = 1e-15)
  expect_equal(r1$TOC$estimate[length(q)], 0, tolerance = 1e-15)
})

test_that("cluster robust rank_average_treatment_effect is consistent", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

  Xc <- rbind(X, X, X, X, X)
  Wc <- c(W, W, W, W, W)
  Yc <- c(Y, Y, Y, Y, Y)
  clust <- rep(1:n, 5)

  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)
  cf.clust <- causal_forest(Xc, Yc, Wc, W.hat = 0.5, clusters = clust, num.trees = 250)
  prio <- runif(n)
  prio.clust <- rep(prio, 5)

  autoc <- rank_average_treatment_effect(cf, prio, target = "AUTOC")
  qini <- rank_average_treatment_effect(cf, prio, target = "QINI")

  autoc.clust <- rank_average_treatment_effect(cf.clust, prio.clust, target = "AUTOC")
  qini.clust <- rank_average_treatment_effect(cf.clust, prio.clust, target = "QINI")

  expect_equal(autoc[["estimate"]], autoc.clust[["estimate"]], tolerance = 0.05)
  expect_equal(autoc[["std.err"]], autoc.clust[["std.err"]], tolerance = 0.02)

  expect_equal(qini[["estimate"]], qini.clust[["estimate"]], tolerance = 0.05)
  expect_equal(qini[["std.err"]], qini.clust[["std.err"]], tolerance = 0.02)
})

test_that("internal bootstrap function `boot_grf` works as expected", {
  n <- 50
  mu <- 1 + rnorm(n)
  clust <- 1:n

  statistic <- function(data, indices) {
    mean(data[indices, 1])
  }

  lm <- lm(mu ~ 1)
  lmt <- lmtest::coeftest(lm, vcov = sandwich::vcovCL, type = "HC", cluster = clust)
  boot.output <- boot_grf(data.frame(mu), statistic, R = 200, clusters = clust, half.sample = FALSE)
  expect_equal(sd(boot.output$t), lmt[1, "Std. Error"], tolerance = 0.02)

  mu.cl <- c(mu, mu, mu, mu)
  clust.cl <- rep(1:n, 4)

  lm.cl <- lm(mu.cl ~ 1)
  lmt.cl <- lmtest::coeftest(lm.cl, vcov = sandwich::vcovCL, type = "HC", cluster = clust.cl)
  boot.output.cl <- boot_grf(data.frame(mu.cl), statistic, R = 200, clusters = clust.cl, half.sample = FALSE)
  expect_equal(sd(boot.output.cl$t), lmt.cl[1, "Std. Error"], tolerance = 0.02)
})

test_that("rank_average_treatment_effect is internally consistent", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  cf <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0.5, num.trees = 250)
  prio1 <- tau
  prio2 <- runif(n)

  set.seed(42)
  rate1 <- rank_average_treatment_effect(cf, cbind(prio1), R = 50)
  set.seed(42)
  rate12 <- rank_average_treatment_effect(cf, cbind(prio1, prio2), R = 50)
  set.seed(42)
  rate2 <- rank_average_treatment_effect(cf, cbind(prio2), R = 50)
  set.seed(42)
  rate21 <- rank_average_treatment_effect(cf, cbind(prio2, prio1), R = 50)
  set.seed(42)
  rate22 <- rank_average_treatment_effect(cf, cbind(prio2, prio2), R = 50)

  expect_equal(rate12$estimate[[1]], rate1$estimate)
  expect_equal(rate12$std.err[[1]], rate1$std.err)
  expect_equal(rate12$estimate[[2]], rate2$estimate)
  expect_equal(rate12$std.err[[2]], rate2$std.err)
  expect_equal(rate12$estimate[[3]], rate1$estimate - rate2$estimate)
  expect_equal(rate12$estimate[[3]], -rate21$estimate[[3]])
  expect_equal(rate12$std.err[[3]], rate21$std.err[[3]])

  expect_equal(rate21$estimate[[1]], rate2$estimate)
  expect_equal(rate21$std.err[[1]], rate2$std.err)
  expect_equal(rate21$estimate[[2]], rate1$estimate)
  expect_equal(rate21$std.err[[2]], rate1$std.err)

  expect_equal(rate22$estimate[[1]], rate22$estimate[[2]])
  expect_equal(rate22$std.err[[1]], rate22$std.err[[2]])
  expect_equal(rate22$estimate[[3]], 0)
  expect_equal(rate22$std.err[[3]], 0)

  priority <- unique(rate12$TOC$priority)
  priority21 <- unique(rate21$TOC$priority)
  expect_equal(rate12$TOC[rate12$TOC$priority == priority[1], ], rate1$TOC)
  expect_equal(as.matrix(rate12$TOC[rate12$TOC$priority == priority[2], -4], rownames.force = FALSE),
               as.matrix(rate2$TOC[, -4]))
  expect_equal(rate12$TOC[rate12$TOC$priority == priority[3], "estimate"],
               rate1$TOC$estimate - rate2$TOC$estimate)
  expect_equal(rate12$TOC[rate12$TOC$priority == priority[3], "estimate"],
               -rate21$TOC[rate21$TOC$priority == priority21[3], "estimate"])
  expect_equal(rate12$TOC[rate12$TOC$priority == priority[3], "std.err"],
               rate21$TOC[rate21$TOC$priority == priority21[3], "std.err"])
})

test_that("rank_average_treatment_effect is internally consistent wrt. subsetting", {
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  wts <- sample(c(1, 2), n, TRUE)
  wts[110:121] <- 0 # This reduces the effective subset to units with wts > 0.
  cf <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0.5, num.trees = 250)
  prio <- runif(n)
  debiasing.wts <- runif(n)

  set.seed(42)
  rate1 <- rank_average_treatment_effect(cf, prio, R = 25, subset = 100:150)
  set.seed(42)
  rate2 <- rank_average_treatment_effect(cf, prio[100:150], R = 25, subset = 100:150)
  set.seed(42)
  rate3 <- rank_average_treatment_effect(cf, prio, R = 25, subset = 100:150, debiasing.weights = debiasing.wts)
  set.seed(42)
  rate4 <- rank_average_treatment_effect(cf, prio, R = 25, subset = 100:150, debiasing.weights = debiasing.wts[100:150])

  expect_equal(rate1, rate2)
  expect_equal(rate3, rate4)
})

test_that("rank_average_treatment_effect has not changed", {
  # Lock in current behavior. A user will expect a given R set.seed to produce the same
  # standard errors. If a future change for example changes how random samples are drawn
  # (a different algo is used, or it is done in parallel) standard errors will change.
  # Even though they are equally valid it breaks user expectations.
  X <- structure(c(3.2, 1, 0.3, -0.7, -1.7, 1.6, 0.3, -0.8, -0.3, -0.5,
                   -1.4, 1.2, -0.8, 0, 0.1, -1.1, -1.6, -0.7, -0.3, 0.3, 0.1, -1.4,
                   -1.8, -1.4, 0.7), .Dim = c(25L, 1L))
  Y <- c(7.2, 4.3, -1.1, -1.9, -3.7, 1.7, 0.4, -0.8, 0.2, -2.4, -2.6,
         2.7, -1.1, 0.3, 0.5, -2.4, -2.9, 0.2, -2.1, -0.1, -0.1, -4, -5.9,
         -3.3, 1.8)
  W <- c(1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L,
         0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 1L)

  cf <- causal_forest(X, Y, W, Y.hat = 0, W.hat = 0.5, num.trees = 250, seed = 42, num.threads = 1)

  set.seed(42)
  rate <- rank_average_treatment_effect(cf, W)

  expect_equal(rate$estimate, -0.84766553531391)
  expect_equal(rate$TOC$estimate, c(-1.05976470588235, -1.05976470588235, -1.05976470588235, -1.05976470588235,
                                    -1.05976470588235, -1.05976470588235, -0.965142857142857, -0.563,
                                    -0.250222222222222, 0))

  expect_equal(rate$std.err, 0.571871384204182)
  expect_equal(rate$TOC$std.err, c(0.795638422842039, 0.795638422842039, 0.795638422842039, 0.787412259502468,
                                   0.755534565676904, 0.714814376377961, 0.573841753475629, 0.394040294086509,
                                   0.188648018385444, 0))
})

test_that("rank_average_treatment_effect.fit works as expected", {
    n <- 500
    p <- 5
    X <- matrix(rnorm(n * p), n, p)
    W <- rbinom(n, 1, 0.5)
    tau <- pmax(X[, 1], 0)
    Y <- tau * W + pmax(X[, 3], 0) + runif(n)
    wts <- sample(1:2, n, TRUE)
    clust <- sample(1:50, n, TRUE)
    cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250, sample.weights = wts, clusters = clust)
    rand <- runif(n)

    set.seed(42)
    rate <- rank_average_treatment_effect(cf, tau)
    set.seed(42)
    rate.fit <- rank_average_treatment_effect.fit(get_scores(cf), tau, sample.weights = wts, clusters = clust)
    expect_equal(rate, rate.fit)

    set.seed(42)
    rate <- rank_average_treatment_effect(cf, list(tau, rand))
    set.seed(42)
    rate.fit <- rank_average_treatment_effect.fit(get_scores(cf), list(tau, rand), sample.weights = wts, clusters = clust)
    expect_equal(rate, rate.fit)
})
