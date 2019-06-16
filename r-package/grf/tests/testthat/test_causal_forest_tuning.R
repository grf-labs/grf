test_that("causal forest tuning decreases prediction error", {
    p = 4
    n = 1000

    X = matrix(2 * runif(n * p) - 1, n, p)
    W = rbinom(n, 1, 0.5)
    TAU = 0.1 * (X[,1] > 0)
    Y = TAU * (W  - 1/2) + 2 * rnorm(n)

    forest = causal_forest(seed=1000, X, Y, W, num.trees = 400, min.node.size = 1, tune.parameters = FALSE)
    preds = predict(forest)
    error = mean((preds$predictions - TAU)^2)

    tuned.forest = causal_forest(seed=1000, X, Y, W, num.trees = 400, tune.parameters = TRUE)
    tuned.preds = predict(tuned.forest)
    tuned.error = mean((tuned.preds$predictions - TAU)^2)

    expect_true(tuned.error < error * 0.75)
})

test_that("causal forest tuning only cross-validates null parameters", {
    p = 6
    n = 100

    X = matrix(2 * runif(n * p) - 1, n, p)
    W = rbinom(n, 1, 0.5)
    TAU = 2 * (X[,1] > 0)
    Y = TAU * (W  - 1/2) + 2 * rnorm(n)

    min.node.size = 5
    imbalance.penalty = 0.42

    tune.output = tune_causal_forest(seed=1000, X, Y, W, 0, 0.5,
                                     min.node.size = min.node.size,
                                     imbalance.penalty = imbalance.penalty)
    tunable.params = tune.output$params

    expect_equal(as.numeric(tunable.params["min.node.size"]), min.node.size)
    expect_equal(as.numeric(tunable.params["imbalance.penalty"]), imbalance.penalty)
})

test_that("local linear causal forest tuning returns lambda and decreases error", {
   p = 6
   n = 1000

   X = matrix(rnorm(n*p), n, p)
   W = rbinom(n, 1, 0.5)
   TAU = 2 * (X[,1] > 0)
   Y = W * TAU + rnorm(n)

   forest = causal_forest(seed=1000, X, Y, W, num.trees = 400)
   tuning.results = tune_ll_causal_forest( forest)
   lambda = tuning.results$lambda.min

   expect_true(is.numeric(lambda))
   expect_true(length(lambda) == 1)

   preds.tuned = predict(forest, linear.correction.variables = 1:p, ll.lambda = lambda)$predictions
   error.tuned = mean((preds.tuned - TAU)^2)

   preds.untuned = predict(forest, linear.correction.variables = 1:p, ll.lambda = 0)$predictions
   error.untuned = mean((preds.untuned - TAU)^2)

   expect_true(error.tuned < 0.75 * error.untuned)
})

test_that("output of tune local linear causal forest is consistent with prediction output", {
  p = 4
  n = 200

  X = matrix(rnorm(n*p), n, p)
  W = rbinom(n, 1, 0.5)
  TAU = 2 * (X[,1] > 0)
  Y = W * TAU + rnorm(n)

  forest = causal_forest(seed=1000, X, Y, W, num.trees = 400)

  tuning.results = tune_ll_causal_forest( forest)

  ll.min = tuning.results$lambdas[1]
  pred.ll.min = predict(forest, linear.correction.variables = 1:p, ll.lambda = ll.min)$predictions
  expect_true(max(abs(tuning.results$oob.predictions[,1] - pred.ll.min)) < 10^-6)

  ll.max = tuning.results$lambdas[length(tuning.results$lambdas)]
  pred.ll.max = predict(forest, linear.correction.variables = 1:p, ll.lambda = ll.max)$predictions
  expect_true(max(abs(tuning.results$oob.predictions[,length(tuning.results$lambdas)] - pred.ll.max)) < 10^-6)
  })
