library(grf)

set.seed(1234)

test_that("calibration is reasonable", {
    n = 400; p = 4
    X = matrix(rnorm(n*p), n, p)
    W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
    Y = pmax(X[,1], 0) * (W - 0.75) + rnorm(n)
    cf = causal_forest(X, Y, W)
    test_calibration(cf)
    expect_true(abs(test_calibration(cf)[1] - 1) <= 0.4)
    
    regf = regression_forest(X, Y)
    test_calibration(regf)
    expect_true(abs(test_calibration(regf)[1]/test_calibration(regf)[2]) <= 3.5)
})
