library(grf)

set.seed(1000)

test_that("average effects are translation invariant", {
	p = 6
	n = 200
	X = matrix(2 * runif(n * p) - 1, n, p)
	W = rbinom(n, 1, 0.5)
	Y = (X[,1] > 0) * (2 * W  - 1) + 2 * rnorm(n)
	Y.plus.1 = Y + 1
	forest.causal = causal_forest(X, Y, W, num.trees = 100, ci.group.size = 4)
	forest.causal.plus.1 = forest.causal
	forest.causal.plus.1$Y.orig = forest.causal$Y.orig + 1
	forest.causal.plus.1$Y.hat = forest.causal$Y.hat + 1

	cate.aipw = average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
	cate.plus.1.aipw = average_treatment_effect(forest.causal.plus.1, target.sample = "all", method = "AIPW")
	expect_equal(cate.aipw, cate.plus.1.aipw)

	cate.tmle = average_treatment_effect(forest.causal, target.sample = "all", method = "TMLE")
	cate.plus.1.tmle = average_treatment_effect(forest.causal.plus.1, target.sample = "all", method = "TMLE")
	expect_equal(cate.tmle, cate.plus.1.tmle)

	catt.aipw = average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
	catt.plus.1.aipw = average_treatment_effect(forest.causal.plus.1, target.sample = "treated", method = "AIPW")
	expect_equal(catt.aipw, catt.plus.1.aipw)

	catt.tmle = average_treatment_effect(forest.causal, target.sample = "treated", method = "TMLE")
	catt.plus.1.tmle = average_treatment_effect(forest.causal.plus.1, target.sample = "treated", method = "TMLE")
	expect_equal(catt.tmle, catt.plus.1.tmle)

	catc.aipw = average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
	catc.plus.1.aipw = average_treatment_effect(forest.causal.plus.1, target.sample = "control", method = "AIPW")
	expect_equal(catc.aipw, catc.plus.1.aipw)

	catc.tmle = average_treatment_effect(forest.causal, target.sample = "control", method = "TMLE")
	catc.plus.1.tmle = average_treatment_effect(forest.causal.plus.1, target.sample = "control", method = "TMLE")
	expect_equal(catc.tmle, catc.plus.1.tmle)

	cape = average_partial_effect(forest.causal)
	cape.plus.1 = average_partial_effect(forest.causal.plus.1)
	expect_true(abs(cape[1] - cape.plus.1[1]) <= 0.005)

	wate = average_treatment_effect(forest.causal, target.sample = "overlap")
	wate.plus.1 = average_treatment_effect(forest.causal.plus.1, target.sample = "overlap")
	expect_true(abs(wate[1] - wate.plus.1[1]) <= 0.005)
})

test_that("average effect estimates are reasonable", {
  p = 6
  n = 1000

  X = matrix(2 * runif(n * p) - 1, n, p)
  eX = 0.25 + 0.5 * (X[,1] > 0)
  W = rbinom(n, 1, eX)
  TAU = 4 * (X[,1] > 0)
  Y =  TAU * (W  - 0.5) + rnorm(n)

  forest.causal = causal_forest(X, Y, W, num.trees = 500, ci.group.size = 1, precompute.nuisance = TRUE)

  cate.aipw = average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  expect_true(abs(cate.aipw[1] - mean(TAU)) <= 0.2)
  expect_true(abs(cate.aipw[1] - mean(TAU)) <= 3 * cate.aipw[2])

  cate.tmle = average_treatment_effect(forest.causal, target.sample = "all", method = "TMLE")
  expect_true(abs(cate.tmle[1] - mean(TAU)) <= 0.2)
  expect_true(abs(cate.tmle[1] - mean(TAU)) <= 3 * cate.tmle[2])

  expect_true(abs(cate.aipw[1] - cate.tmle[1]) <= 0.01)
  expect_true(abs(cate.aipw[2] - cate.tmle[2]) <= 0.01)

  catt.aipw = average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  expect_true(abs(catt.aipw[1] - mean(TAU[W==1])) <= 0.2)
  expect_true(abs(catt.aipw[1] - mean(TAU[W==1])) <= 3 * catt.aipw[2])

  catt.tmle = average_treatment_effect(forest.causal, target.sample = "treated", method = "TMLE")
  expect_true(abs(catt.tmle[1] - mean(TAU[W==1])) <= 0.2)
  expect_true(abs(catt.tmle[1] - mean(TAU[W==1])) <= 3 * catt.tmle[2])

  expect_true(abs(catt.aipw[1] - catt.tmle[1]) <= 0.05)
  expect_true(abs(catt.aipw[2] - catt.tmle[2]) <= 0.05)

  catc.aipw = average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  expect_true(abs(catc.aipw[1] - mean(TAU[W==0])) <= 0.2)
  expect_true(abs(catc.aipw[1] - mean(TAU[W==0])) <= 3 * catc.aipw[2])

  catc.tmle = average_treatment_effect(forest.causal, target.sample = "control", method = "TMLE")
  expect_true(abs(catc.tmle[1] - mean(TAU[W==0])) <= 0.2)
  expect_true(abs(catc.tmle[1] - mean(TAU[W==0])) <= 3 * catc.tmle[2])

  expect_true(abs(catc.aipw[1] - catc.tmle[1]) <= 0.05)
  expect_true(abs(catc.aipw[2] - catc.tmle[2]) <= 0.05)

  cape = average_partial_effect(forest.causal)
  expect_true(abs(cape[1] - mean(TAU)) <= 0.2)
  expect_true(abs(cape[1] - mean(TAU)) <= 3 * cape[2])

  expect_true(abs(cate.aipw[1] - cape[1]) <= 0.05)
  expect_true(abs(cate.aipw[2] - cape[2]) <= 0.05)

  wate = average_treatment_effect(forest.causal, target.sample = "overlap")
  tau.overlap = sum(eX * (1 - eX) * TAU) / sum(eX * (1 - eX))
  expect_true(abs(wate[1] - tau.overlap) <= 0.2)
  expect_true(abs(wate[1] - tau.overlap) <= 3 * wate[2])
})

test_that("average treatment effects larger example works", {

  n = 4000
  p = 10

  X = matrix(runif(n * p), n, p)
  E = (0.4 + dbeta(X[,2], 2, 4)) / 4
  W = rbinom(n, 1, E)
  M = 2 * X[,2] - 1
  TAU = (1 + 1/(1 + exp(-20 * (X[,1] - 0.3)))) * (1 + 1/(1 + exp(-20 * (X[,2] - 0.3))))
  Y = M + (W - 0.5) * TAU + rnorm(n)

  forest.causal = causal_forest(X, Y, W, num.trees = 1000, ci.group.size = 1, precompute.nuisance = TRUE)

  cate.aipw = average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  expect_true(abs(cate.aipw[1] - mean(TAU)) <= 3 * cate.aipw[2])

  cate.tmle = average_treatment_effect(forest.causal, target.sample = "all", method = "TMLE")
  expect_true(abs(cate.tmle[1] - mean(TAU)) <= 3 * cate.tmle[2])

  catt.aipw = average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  expect_true(abs(catt.aipw[1] - mean(TAU[W==1])) <= 3 * catt.aipw[2])

  catt.tmle = average_treatment_effect(forest.causal, target.sample = "treated", method = "TMLE")
  expect_true(abs(catt.tmle[1] - mean(TAU[W==1])) <= 3 * catt.tmle[2])

  catc.aipw = average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  expect_true(abs(catc.aipw[1] - mean(TAU[W==0])) <= 3 * catc.aipw[2])

  catc.tmle = average_treatment_effect(forest.causal, target.sample = "control", method = "TMLE")
  expect_true(abs(catc.tmle[1] - mean(TAU[W==0])) <= 3 * catc.tmle[2])
})

test_that("average partial effects larger example works", {

  n = 4000
  p = 10

  X = matrix(runif(n * p), n, p)
  E = (0.4 + dbeta(X[,2], 2, 4)) / 4
  W = rbinom(n, 1, E) + 0.2 * rnorm(n)
  M = 2 * X[,2] - 1
  TAU = (1 + 1/(1 + exp(-20 * (X[,1] - 0.3)))) * (1 + 1/(1 + exp(-20 * (X[,2] - 0.3))))
  Y = M + (W - 0.5) * TAU + rnorm(n)

  forest.causal = causal_forest(X, Y, W, num.trees = 1000, ci.group.size = 1, precompute.nuisance = TRUE)

  cape = average_partial_effect(forest.causal)
  expect_true(abs(cape[1] - mean(TAU)) <= 0.2)
  expect_true(abs(cape[1] - mean(TAU)) <= 3 * cape[2])

})

test_that("average treatment effect with overlap: larger example works", {

  n = 4000
  p = 10

  X = matrix(rnorm(n * (p)), n, p)
  eX = 1/(1 + exp(-10 * X[,2]))
  W = rbinom(n, 1, eX)
  M = X[,2]
  TAU = (1 + X[,2])^2
  Y = M + (W - 0.5) * TAU + rnorm(n)

  forest.causal = causal_forest(X, Y, W, num.trees = 1000, ci.group.size = 1, precompute.nuisance = TRUE)

  wate = average_treatment_effect(forest.causal, target.sample = "overlap")
  tau.overlap = sum(eX * (1 - eX) * TAU) / sum(eX * (1 - eX))
  expect_true(abs(wate[1] - tau.overlap) <= 0.2)
  expect_true(abs(wate[1] - tau.overlap) <= 3 * wate[2])

})

test_that("cluster robust average effects are consistent", {
  p = 6
  n = 400

  X = matrix(2 * runif(n * p) - 1, n, p)
  W = rbinom(n, 1, 0.5)
  Y = (X[,1] > 0) * (2 * W  - 1) + 2 * rnorm(n)

  Xc = rbind(X, X, X, X)
  Wc = c(W, W, W, W)
  Yc = c(Y, Y, Y, Y)
  clust = c(1:n, 1:n, 1:n, 1:n)

  forest.causal = causal_forest(X, Y, W, num.trees = 4000, ci.group.size = 4)
  forest.causal.clust = causal_forest(Xc, Yc, Wc, num.trees = 4000, ci.group.size = 4, clusters = clust)

  cate.aipw = average_treatment_effect(forest.causal, target.sample = "all", method = "AIPW")
  cate.clust.aipw = average_treatment_effect(forest.causal.clust, target.sample = "all", method = "AIPW")
  expect_true(abs(cate.aipw[1] - cate.clust.aipw[1]) <= 0.05)
  expect_true(abs(cate.aipw[2] - cate.clust.aipw[2]) <= 0.005)

  cate.tmle = average_treatment_effect(forest.causal, target.sample = "all", method = "TMLE")
  cate.clust.tmle = average_treatment_effect(forest.causal.clust, target.sample = "all", method = "TMLE")
  expect_true(abs(cate.tmle[1] - cate.clust.tmle[1]) <= 0.05)
  expect_true(abs(cate.tmle[2] - cate.clust.tmle[2]) <= 0.005)

  catt.aipw = average_treatment_effect(forest.causal, target.sample = "treated", method = "AIPW")
  catt.clust.aipw = average_treatment_effect(forest.causal.clust, target.sample = "treated", method = "AIPW")
  expect_true(abs(catt.aipw[1] - catt.clust.aipw[1]) <= 0.1)
  expect_true(abs(catt.aipw[2] - catt.clust.aipw[2]) <= 0.005)

  catt.tmle = average_treatment_effect(forest.causal, target.sample = "treated", method = "TMLE")
  catt.clust.tmle = average_treatment_effect(forest.causal.clust, target.sample = "treated", method = "TMLE")
  expect_true(abs(catt.tmle[1] - catt.clust.tmle[1]) <= 0.1)
  expect_true(abs(catt.tmle[2] - catt.clust.tmle[2]) <= 0.005)

  catc.aipw = average_treatment_effect(forest.causal, target.sample = "control", method = "AIPW")
  catc.clust.aipw = average_treatment_effect(forest.causal.clust, target.sample = "control", method = "AIPW")
  expect_true(abs(catc.aipw[1] - catc.clust.aipw[1]) <= 0.1)
  expect_true(abs(catc.aipw[2] - catc.clust.aipw[2]) <= 0.005)

  catc.tmle = average_treatment_effect(forest.causal, target.sample = "control", method = "TMLE")
  catc.clust.tmle = average_treatment_effect(forest.causal.clust, target.sample = "control", method = "TMLE")
  expect_true(abs(catc.tmle[1] - catc.clust.tmle[1]) <= 0.1)
  expect_true(abs(catc.tmle[2] - catc.clust.tmle[2]) <= 0.005)

  cape = average_partial_effect(forest.causal)
  cape.clust = average_partial_effect(forest.causal.clust)
  expect_true(abs(cape[1] - cape.clust[1]) <= 0.1)
  expect_true(abs(cape[2] - cape.clust[2]) <= 0.01)

  wate = average_treatment_effect(forest.causal, target.sample = "overlap")
  wate.clust = average_treatment_effect(forest.causal.clust, target.sample = "overlap")
  expect_true(abs(wate[1] - wate.clust[1]) <= 0.05)
  expect_true(abs(wate[2] - wate.clust[2]) <= 0.005)
})