# *** Comparison methods ***

#' Compute E[T | X]
#'
#' @param S.hat The estimated survival curve.
#' @param Y.grid The time values corresponding to S.hat.
#' @return A vector of expected values.
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))

  c(cbind(1, S.hat) %*% grid.diff)
}

# "SRC1"
estimate_rfsrc_X_W = function(data, data.test) {
  df = data.frame(Y=data$Y, D=data$D, cbind(x=data$X, w=data$W))
  ntree = 500
  fit = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df)
  n.test = nrow(data.test$X)
  W1 = rep(1, n.test)
  W0 = rep(0, n.test)
  p1 = predict(fit, data.frame(cbind(data.test$X, w=W1)))
  p0 = predict(fit, data.frame(cbind(data.test$X, w=W0)))

  expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
}

# "SRC2"
estimate_rfsrc_XW_W = function(data, data.test) {
  df = data.frame(Y=data$Y, D=data$D, cbind(x=data$X, xw=data$X*data$W, w=data$W))
  ntree = 500
  fit = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df)
  n.test = nrow(data.test$X)
  W1 = rep(1, n.test)
  W0 = rep(0, n.test)
  p1 = predict(fit, data.frame(cbind(data.test$X, data.test$X*W1, w=W1)))
  p0 = predict(fit, data.frame(cbind(data.test$X, data.test$X*W0, w=W0)))

  expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
}

# "VT"
estimate_rfsrc_twin = function(data, data.test) {
  ntree = 500
  W = data$W
  df1 = data.frame(Y=data$Y[W==1], D=data$D[W==1], data$X[W==1, ])
  df0 = data.frame(Y=data$Y[W==0], D=data$D[W==0], data$X[W==0, ])
  fit1 = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df1)
  fit0 = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df0)
  p1 = predict(fit1, data.frame(data.test$X))
  p0 = predict(fit0, data.frame(data.test$X))

  expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
}

# "CSF"
estimate_grf = function(data, data.test) {
  fit = causal_survival_forest(data$X, data$Y, data$W, data$D)
  p = predict(fit, data.test$X)

  p$predictions
}

# "IPCW"
estimate_IPCW_grf = function(data, data.test) {
  sf.censor = survival_forest(data$X, data$Y, 1 - data$D, prediction.type = "Nelson-Aalen",
                              num.trees = 500)
  C.hat = predict(sf.censor)$predictions
  Y.relabeled = sf.censor$Y.relabeled
  C.Y.hat = rep(1, nrow(data$X)) # (for events before the first failure, C.Y.hat is one)
  C.Y.hat[Y.relabeled != 0] = C.hat[cbind(1:nrow(data$X), Y.relabeled)]
  sample.weights = 1 / C.Y.hat
  subset = data$D == 1
  cf = causal_forest(data$X[subset, ], data$Y[subset], data$W[subset], sample.weights = sample.weights[subset])
  p = predict(cf, data.test$X)

  p$predictions
}
