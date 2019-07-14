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
	
	cate.aipw = estimate_average_effect(forest.causal, target.sample = "all", method = "AIPW")
	cate.plus.1.aipw = estimate_average_effect(forest.causal.plus.1, target.sample = "all", method = "AIPW")
	expect_equal(cate.aipw, cate.plus.1.aipw)
	
	cate.tmle = estimate_average_effect(forest.causal, target.sample = "all", method = "TMLE")
	cate.plus.1.tmle = estimate_average_effect(forest.causal.plus.1, target.sample = "all", method = "TMLE")
	expect_equal(cate.tmle, cate.plus.1.tmle)
	
	catt.aipw = estimate_average_effect(forest.causal, target.sample = "treated", method = "AIPW")
	catt.plus.1.aipw = estimate_average_effect(forest.causal.plus.1, target.sample = "treated", method = "AIPW")
	expect_equal(catt.aipw, catt.plus.1.aipw)
	
	catt.tmle = estimate_average_effect(forest.causal, target.sample = "treated", method = "TMLE")
	catt.plus.1.tmle = estimate_average_effect(forest.causal.plus.1, target.sample = "treated", method = "TMLE")
	expect_equal(catt.tmle, catt.plus.1.tmle)
	
	catc.aipw = estimate_average_effect(forest.causal, target.sample = "control", method = "AIPW")
	catc.plus.1.aipw = estimate_average_effect(forest.causal.plus.1, target.sample = "control", method = "AIPW")
	expect_equal(catc.aipw, catc.plus.1.aipw)
	
	catc.tmle = estimate_average_effect(forest.causal, target.sample = "control", method = "TMLE")
	catc.plus.1.tmle = estimate_average_effect(forest.causal.plus.1, target.sample = "control", method = "TMLE")
	expect_equal(catc.tmle, catc.plus.1.tmle)
})

test_that("average effect estimates are reasonable", {
  p = 6
  n = 1000
  
  X = matrix(2 * runif(n * p) - 1, n, p)
  W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
  TAU = 4 * (X[,1] > 0)
  Y =  TAU * (W  - 0.5) + rnorm(n)
  
  forest.causal = causal_forest(X, Y, W, num.trees = 500, ci.group.size = 1, precompute.nuisance = TRUE)
  
  cate.aipw = estimate_average_effect(forest.causal, target.sample = "all", method = "AIPW")
  expect_true(abs(cate.aipw[1] - mean(TAU)) <= 0.2)
  expect_true(abs(cate.aipw[1] - mean(TAU)) <= 3 * cate.aipw[2])
  
  cate.tmle = estimate_average_effect(forest.causal, target.sample = "all", method = "TMLE")
  expect_true(abs(cate.tmle[1] - mean(TAU)) <= 0.2)
  expect_true(abs(cate.tmle[1] - mean(TAU)) <= 3 * cate.tmle[2])
  
  catt.aipw = estimate_average_effect(forest.causal, target.sample = "treated", method = "AIPW")
  expect_true(abs(catt.aipw[1] - mean(TAU[W==1])) <= 0.2)
  expect_true(abs(catt.aipw[1] - mean(TAU[W==1])) <= 3 * catt.aipw[2])
  
  catt.tmle = estimate_average_effect(forest.causal, target.sample = "treated", method = "TMLE")
  expect_true(abs(catt.tmle[1] - mean(TAU[W==1])) <= 0.2)
  expect_true(abs(catt.tmle[1] - mean(TAU[W==1])) <= 3 * catt.tmle[2])
  
  catc.aipw = estimate_average_effect(forest.causal, target.sample = "control", method = "AIPW")
  expect_true(abs(catc.aipw[1] - mean(TAU[W==0])) <= 0.2)
  expect_true(abs(catc.aipw[1] - mean(TAU[W==0])) <= 3 * catc.aipw[2])
  
  catc.tmle = estimate_average_effect(forest.causal, target.sample = "control", method = "TMLE")
  expect_true(abs(catc.tmle[1] - mean(TAU[W==0])) <= 0.2)
  expect_true(abs(catc.tmle[1] - mean(TAU[W==0])) <= 3 * catc.tmle[2])
})

test_that("average effects larger example works", {
  
  n = 4000
  p = 10
  
  X = matrix(runif(n * p), n, p)
  E = (0.4 + dbeta(X[,2], 2, 4)) / 4
  W = rbinom(n, 1, E)
  M = 2 * X[,2] - 1
  TAU = (1 + 1/(1 + exp(-20 * (X[,1] - 0.3)))) * (1 + 1/(1 + exp(-20 * (X[,2] - 0.3))))
  Y = M + (W - 0.5) * TAU + rnorm(n)
  
  forest.causal = causal_forest(X, Y, W, num.trees = 1000, ci.group.size = 1, precompute.nuisance = TRUE)
  
  cate.aipw = estimate_average_effect(forest.causal, target.sample = "all", method = "AIPW")
  expect_true(abs(cate.aipw[1] - mean(TAU)) <= 3 * cate.aipw[2])
  
  cate.tmle = estimate_average_effect(forest.causal, target.sample = "all", method = "TMLE")
  expect_true(abs(cate.tmle[1] - mean(TAU)) <= 3 * cate.tmle[2])
  
  catt.aipw = estimate_average_effect(forest.causal, target.sample = "treated", method = "AIPW")
  expect_true(abs(catt.aipw[1] - mean(TAU[W==1])) <= 3 * catt.aipw[2])
  
  catt.tmle = estimate_average_effect(forest.causal, target.sample = "treated", method = "TMLE")
  expect_true(abs(catt.tmle[1] - mean(TAU[W==1])) <= 3 * catt.tmle[2])
  
  catc.aipw = estimate_average_effect(forest.causal, target.sample = "control", method = "AIPW")
  expect_true(abs(catc.aipw[1] - mean(TAU[W==0])) <= 3 * catc.aipw[2])
  
  catc.tmle = estimate_average_effect(forest.causal, target.sample = "control", method = "TMLE")
  expect_true(abs(catc.tmle[1] - mean(TAU[W==0])) <= 3 * catc.tmle[2])
})


