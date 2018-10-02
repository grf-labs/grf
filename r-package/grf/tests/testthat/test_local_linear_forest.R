library(grf)

set.seed(1234)

test_that("local linear prediction gives reasonable estimates", {
    f = function(x){x[1] + 2*x[2] + 2*x[3]**2}
    n = 600
    p = 5
    X = matrix(rnorm(n*p), n, p)
    MU = apply(X, FUN=f, MARGIN=1)
    Y = MU + rnorm(n)

    forest = regression_forest(X, Y)
    preds.grf.oob = predict(forest)
    preds.ll.oob = predict(forest, linear.correction.variables = 1:p, ll.lambda = 0)

    mse.grf.oob = mean((preds.grf.oob$predictions - MU)^2)
    mse.ll.oob = mean((preds.ll.oob$predictions - MU)^2)

    expect_true(mse.ll.oob < 1)
    expect_true(mse.ll.oob < mse.grf.oob / 2)

    X.test = matrix(rnorm(n*p), n, p)
    MU.test = apply(X.test, FUN=f, MARGIN=1)

    preds.grf = predict(forest, X.test)
    preds.ll = predict(forest, X.test, linear.correction.variables = 1:p, ll.lambda = 0.1)

    mse.grf = mean((preds.grf$predictions - MU.test)^2)
    mse.ll = mean((preds.ll$predictions - MU.test)^2)

    expect_true(mse.ll < 1)
    expect_true(mse.ll < mse.grf / 1.5)
})

test_that("linear correction variables function as expected", {
    f = function(x){x[1] + 2*x[2] + 2*x[3]**2}
    n = 400
    p = 20
    X = matrix(rnorm(n*p), n, p)
    MU = apply(X, FUN=f, MARGIN=1)
    Y = MU + rnorm(n)

    forest = regression_forest(X, Y)
    preds = predict(forest, linear.correction.variables = 1:20)
    mse = mean((preds$predictions - MU)^2)

    preds.selected = predict(forest, linear.correction.variables = 1:3)
    mse.selected = mean((preds.selected$predictions - MU)^2)

    expect_true(mse.selected < mse / 1.5)
})

test_that("local linear forest tuning always returns lambda and decreases prediction error", {
    n = 1000
    p = 20
    sigma = 5

    mu = function(x){log(1+exp(6*x[1]))}

    X = matrix(runif(n*p,-1,1), nrow = n)
    truth = apply(X, FUN = mu, MARGIN = 1)
    Y = truth + sigma*rnorm(n)

    forest = local_linear_forest(X, Y)

    lambda = tune_local_linear_forest(forest)$lambda.min
    expect_true(is.numeric(lambda))
    expect_true(length(lambda) == 1)

    preds.tuned = predict(forest, linear.correction.variables = 1:20, ll.lambda = lambda)$predictions
    mse.tuned = mean((preds.tuned - truth)^2)

    preds.untuned = predict(forest, linear.correction.variables = 1:20, ll.lambda = 0.1)$predictions
    mse.untuned = mean((preds.untuned - truth)^2)

    expect_true(mse.tuned < 0.75 * mse.untuned)
})

test_that("default local linear forest predict and regression forest predict with local.linear = TRUE are the same", {
    n = 1000
    p = 5
    sigma = 1

    mu = function(x){log(1+exp(6*x[1]))}

    X = matrix(runif(n*p,-1,1), nrow = n)
    truth = apply(X, FUN = mu, MARGIN = 1)
    Y = truth + sigma*rnorm(n)

    forest = regression_forest(X, Y)
    preds = predict(forest, linear.correction.variables = 1:5, lambda = 0.1)$predictions

    ll.forest = local_linear_forest(X, Y)
    ll.preds = predict(ll.forest, ll.lambda = 0.1)$predictions

    average.difference = mean((ll.preds - preds)**2)

    expect_true(average.difference < 0.1)
})

test_that("local linear predict returns local linear predictions even without tuning parameters", {
    n = 200
    p = 5
    sigma = 1

    mu = function(x){log(1+exp(6*x[1]))}

    X = matrix(runif(n*p,-1,1), nrow = n)
    truth = apply(X, FUN = mu, MARGIN = 1)
    Y = truth + sigma*rnorm(n)

    forest = local_linear_forest(X, Y)
    preds = predict(forest)

    ll.indicator = !is.null(preds$ll.lambda)
    expect_true(ll.indicator)
})

test_that("local linear confidence intervals have reasonable coverage", {
    mu = function(x){log(1+exp(6*x))}

    n = 500
    p = 20
    sigma = sqrt(20)

    X = matrix(runif(n*p,-1,1), nrow = n)
    truth = mu(X[,1])
    Y = truth + sigma*rnorm(n)

    forest = local_linear_forest(X, Y, num.trees = 1000, ci.group.size = 5)
    preds = predict(forest, linear.correction.variables = 1, tune.lambda = TRUE, estimate.variance = TRUE)

    expect_true(all(preds$variance.estimates > 0))

    df = data.frame(predictions= preds$predictions,
                    upper = preds$predictions + 1.96*sqrt(preds$variance.estimates),
                    lower = preds$predictions - 1.96*sqrt(preds$variance.estimates))

    percent_llf = mean(sapply(1:n, function(i){
        ifelse(df$lower[i] <= truth[i] && truth[i] <= df$upper[i], 1, 0)
    }))

    expect_true(percent_llf > 0.94)
})