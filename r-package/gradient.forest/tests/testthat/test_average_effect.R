library(gradient.forest)

set.seed(4321)

test_that("average effects are translation invariant", {
	p = 6
	n = 200

	X = matrix(2 * runif(n * p) - 1, n, p)
	W = rbinom(n, 1, 0.5)
	Y = (X[,1] > 0) * (2 * W  - 1) + 2 * rnorm(n)
	Y.plus.1 = Y + 1

	forest.causal = causal.forest(X, Y, W, num.trees = 100, ci.group.size = 4)
	forest.causal.plus.1 = forest.causal
	forest.causal.plus.1$Y.orig = forest.causal$Y.orig + 1
	forest.causal.plus.1$Y.hat = forest.causal$Y.hat + 1
	
	cate = estimate.average.effect(forest.causal, target.sample = "all")
	cate.plus.1 = estimate.average.effect(forest.causal.plus.1, target.sample = "all")
	expect_equal(cate, cate.plus.1)
	
	catt = estimate.average.effect(forest.causal, target.sample = "treated")
	catt.plus.1 = estimate.average.effect(forest.causal.plus.1, target.sample = "treated")
	expect_equal(catt, catt.plus.1)
	
	catc = estimate.average.effect(forest.causal, target.sample = "control")
	catc.plus.1 = estimate.average.effect(forest.causal.plus.1, target.sample = "control")
	expect_equal(catc, catc.plus.1)
})

test_that("average effect estimates are reasonable", {
  p = 6
  n = 1000
  
  X = matrix(2 * runif(n * p) - 1, n, p)
  W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
  TAU = 4 * (X[,1] > 0)
  Y =  TAU * (W  - 0.5) + rnorm(n)
  
  forest.causal = causal.forest(X, Y, W, num.trees = 500, ci.group.size = 1, precompute.nuisance = TRUE)
  
  cate = estimate.average.effect(forest.causal, target.sample = "all")
  expect_true(abs(cate[1] - mean(TAU)) <= 0.2)
  expect_true(abs(cate[1] - mean(TAU)) <= 3 * cate[2])
  
  catt = estimate.average.effect(forest.causal, target.sample = "treated")
  expect_true(abs(catt[1] - mean(TAU[W==1])) <= 0.2)
  expect_true(abs(catt[1] - mean(TAU[W==1])) <= 3 * catt[2])
  
  catc = estimate.average.effect(forest.causal, target.sample = "control")
  expect_true(abs(catc[1] - mean(TAU[W==0])) <= 0.2)
  expect_true(abs(catc[1] - mean(TAU[W==0])) <= 3 * catc[2])
})

test_that("average effects larger example works", {
  
  n = 4000
  p = 10
  
  X = matrix(runif(n * p), n, p)
  E = (1 + dbeta(X[,2], 2, 4)) / 4
  W = rbinom(n, 1, E)
  M = 2 * X[,2] - 1
  TAU = (1 + 1/(1 + exp(-20 * (X[,1] - 0.3)))) * (1 + 1/(1 + exp(-20 * (X[,2] - 0.3))))
  Y = M + (W - 0.5) * TAU + rnorm(n)
  
  forest.causal = causal.forest(X, Y, W, num.trees = 1000, ci.group.size = 1, precompute.nuisance = TRUE)
  
  cate = estimate.average.effect(forest.causal, target.sample = "all")
  expect_true(abs(cate[1] - mean(TAU)) <= 3 * cate[2])
  
  catt = estimate.average.effect(forest.causal, target.sample = "treated")
  expect_true(abs(catt[1] - mean(TAU[W==1])) <= 3 * catt[2])
  
  catc = estimate.average.effect(forest.causal, target.sample = "control")
  expect_true(abs(catc[1] - mean(TAU[W==0])) <= 3 * catc[2])
})


