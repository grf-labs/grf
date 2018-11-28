rm(list = ls())
library(grf)

generate_data = function(setup) {
  if (setup == 1) {
    n = 600
    p = 8
    X = matrix(runif(n * p), n, p)
    propensity = (1 + dbeta(X[,3], 2, 4)) / 4
    tau = (1 + 1/(1 + exp(-20 * (X[,1] - 0.3)))) * (1 + 1/(1 + exp(-20 * (X[,2] - 0.3))))
    W = rbinom(n, 1, 0.9)
    Y = 2 * X[,3] - 1 + (W - 0.5) * tau + rnorm(n)
  } else if (setup == 2) {
    n = 1600
    p = 6
    k = 3
    X = matrix(rnorm(n*p), n, p)
    a = rowMeans(X[,1:k]) * sqrt(k)
    x = sign(a) * a^2
    tau = 1/(1 + exp(-X[,p]))
    propensity = pmax(0.05, pmin(1/(1 + exp(-x)), 0.95))
    mu = x - tau * (propensity - 0.5)
    W = rbinom(n, 1, propensity)
    Y = mu + W * tau + rnorm(n)
  } else if (setup == 3) {
    n = 1000
    p = 10
    X = matrix(runif(n * p), n, p)
    tau = (1 + 1/(1 + exp(-20 * (X[,1] - 0.3)))) * (1 + 1/(1 + exp(-20 * (X[,2] - 0.3))))
    W = rbinom(n, 1, 0.9)
    Y = (W - 0.9) * tau + rnorm(n)
  } else if (setup == 4) {
    n = 1600
    p = 20
    X = matrix(runif(n * p), n, p)
    propensity = (1 + dbeta(X[,3], 2, 4)) / 4
    tau = rep(0, n)
    W = rbinom(n, 1, 0.9)
    Y = 2 * X[,3] - 1 + (W - 0.5) * tau + rnorm(n)
  } else if (setup == 5) {
    n = 4000
    p = 10
    X = matrix(rnorm(n * p), n, p)
    W = rbinom(n, 1, 0.1)
    tau = 0.2 * (X[,3] > 0)
    Y = X[,1] + X[,2] + tau * W + rnorm(n)
  }
  list(X=X, W=W, Y=Y, tau=tau, n=n, p=p)
}

evaluate_method = function(estimate_tau, setup = 1) {
  data = generate_data(setup)
  tau.hat = estimate_tau(data$X, data$Y, data$W)
  plot(data$tau, tau.hat)
  abline(0, 1)
  sqrt(mean((tau.hat - data$tau)^2))
}

make_causal_forest = function(stabilize.splits, min.node.size,
                              alpha, imbalance.penalty) {
  function(X, Y, W) {
    cf = causal_forest(X, Y, W, stabilize.splits = stabilize.splits, alpha = alpha,
                       min.node.size = min.node.size, imbalance.penalty = imbalance.penalty)
    cf.pred = predict(cf)
    print(cf)
    cf.pred$predictions
  }
}

res.untuned.unstab = sapply(1:5, function(setup) {
  evaluate_method(function(X, Y, W) {
    cf = causal_forest(X, Y, W, tune.parameters = FALSE, stabilize.splits = FALSE)
    cf.pred = predict(cf)
    cf.pred$predictions
  }, setup)
})
res.untuned.unstab

res.tuned.unstab = sapply(1:5, function(setup) {
  evaluate_method(function(X, Y, W) {
    cf = causal_forest(X, Y, W, tune.parameters = TRUE, stabilize.splits = FALSE)
    cf.pred = predict(cf)
    cf.pred$predictions
  }, setup)
})
res.tuned.unstab

res.untuned.stab = sapply(1:5, function(setup) {
  evaluate_method(function(X, Y, W) {
    cf = causal_forest(X, Y, W, min.node.size = 5, tune.parameters = FALSE, stabilize.splits = TRUE)
    cf.pred = predict(cf)
    cf.pred$predictions
  }, setup)
})
res.untuned.stab

res.tuned.stab = sapply(1:5, function(setup) {
  evaluate_method(function(X, Y, W) {
    cf = causal_forest(X, Y, W, tune.parameters = TRUE, stabilize.splits = TRUE)
    cf.pred = predict(cf)
    cf.pred$predictions
  }, setup)
})
res.tuned.stab

res.untuned.unstab
res.tuned.unstab
res.untuned.stab
res.tuned.stab


# res.ip = outer(c(0, 0.25, 0.5, 2, 4), 1:5,
#             FUN = Vectorize(function(imbalance.penalty, setup) {
#               evaluate_method(make_causal_forest(
#                 stabilize.splits = TRUE,
#                 min.node.size = 5,
#                 alpha = 0.05,
#                 imbalance.penalty = imbalance.penalty),
#                 setup)
#             }))
# 
# res.mns = outer(c(0, 1, 2, 4, 8, 16, 32), 1:5,
#                 FUN = Vectorize(function(mns, setup) {
#                   evaluate_method(make_causal_forest(
#                     stabilize.splits = TRUE,
#                     min.node.size = mns,
#                     alpha = 0.05,
#                     imbalance.penalty = 0),
#                     setup)
#                 }))
# colnames(res.mns) = sapply(1:5, function(ii) paste("setup", ii))
# rownames(res.mns) = sapply(c(0, 1, 2, 4, 8, 16, 32), function(ii) paste("min. node size", ii))
