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
SRC1 = function(data, data.test) {
  df = data.frame(Y=data$Y, D=data$D, cbind(x=data$X, w=data$W))
  ntree = 500
  fit = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df)
  n.test = nrow(data.test$X)
  W1 = rep(1, n.test)
  W0 = rep(0, n.test)
  p1 = predict(fit, data.frame(cbind(data.test$X, w=W1)))
  p0 = predict(fit, data.frame(cbind(data.test$X, w=W0)))
  horizonS.index <- findInterval(data$y0, p1$time.interest)

  pp = expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
  pp.prob = p1$survival[, horizonS.index] - p0$survival[, horizonS.index]

  list(pp = pp, pp.prob = pp.prob)
}

# "SRC2"
SRC2 = function(data, data.test) {
  df = data.frame(Y=data$Y, D=data$D, cbind(x=data$X, xw=data$X*data$W, w=data$W))
  ntree = 500
  fit = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df)
  n.test = nrow(data.test$X)
  W1 = rep(1, n.test)
  W0 = rep(0, n.test)
  p1 = predict(fit, data.frame(cbind(data.test$X, data.test$X*W1, w=W1)))
  p0 = predict(fit, data.frame(cbind(data.test$X, data.test$X*W0, w=W0)))
  horizonS.index <- findInterval(data$y0, p1$time.interest)

  pp = expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
  pp.prob = p1$survival[, horizonS.index] - p0$survival[, horizonS.index]

  list(pp = pp, pp.prob = pp.prob)
}

# "VT"
VT = function(data, data.test) {
  ntree = 500
  W = data$W
  df1 = data.frame(Y=data$Y[W==1], D=data$D[W==1], data$X[W==1, ])
  df0 = data.frame(Y=data$Y[W==0], D=data$D[W==0], data$X[W==0, ])
  fit1 = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df1)
  fit0 = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df0)
  p1 = predict(fit1, data.frame(data.test$X))
  p0 = predict(fit0, data.frame(data.test$X))
  horizonS.index1 <- findInterval(data$y0, p1$time.interest)
  horizonS.index0 <- findInterval(data$y0, p0$time.interest)

  pp = expected_survival(p1$survival, p1$time.interest) -
    expected_survival(p0$survival, p0$time.interest)
  pp.prob = p1$survival[, horizonS.index1] - p0$survival[, horizonS.index0]

  list(pp = pp, pp.prob = pp.prob)
}

# "CSF"
CSF = function(data, data.test) {
  fit = causal_survival_forest(data$X, data$Y, data$W, data$D, horizon = data$Y.max)
  fit.prob = causal_survival_forest(data$X, data$Y, data$W, data$D,
                                    target = "survival.probability", horizon = data$y0)
  p = predict(fit, data.test$X)
  p.prob = predict(fit.prob, data.test$X)

  list(pp = p$predictions, pp.prob = p.prob$predictions)
}

# "IPCW"
IPCW = function(data, data.test) {
  sf.censor = survival_forest(cbind(data$X, data$W), data$Y, 1 - data$D,
                              prediction.type = "Nelson-Aalen", num.trees = 500)
  C.hat = predict(sf.censor)$predictions
  Y.relabeled = sf.censor$Y.relabeled
  C.Y.hat = rep(1, nrow(data$X)) # (for events before the first failure, C.Y.hat is one)
  # Pick out S_C(Yi, X) from the estimated survival curve
  C.Y.hat[Y.relabeled != 0] = C.hat[cbind(1:nrow(data$X), Y.relabeled)]
  sample.weights = 1 / C.Y.hat
  subset = data$D == 1
  cf = causal_forest(data$X[subset, ], data$Y[subset], data$W[subset], sample.weights = sample.weights[subset])

  subset = data$D == 1 | data$Y > data$y0
  horizonC.index = findInterval(data$y0, sf.censor$failure.times)
  if (horizonC.index != 0) {
    C.Y.hat[data$Y > data$y0] = C.hat[data$Y > data$y0, horizonC.index]
  }
  sample.weights.pr = 1 / C.Y.hat
  cf.prob = causal_forest(data$X[subset, ], as.numeric(data$Y[subset] > data$y0), data$W[subset], sample.weights = sample.weights.pr[subset])
  p = predict(cf, data.test$X)
  p.prob = predict(cf.prob, data.test$X)

  list(pp = p$predictions, pp.prob = p.prob$predictions)
}
