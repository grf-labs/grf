# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
rm(list = ls())
set.seed(42)
library(survival)
library(randomForestSRC) # 2.9.3
library(ranger) # 0.12.1
library(grf) # 1.2.0.0
source('generate_survival_data.R')

#' Compute E[T | X]
#'
#' @param S.hat The estimated survival curve, a (num.samples * num.events) matrix
#' @param Y.grid The num.events time values corresponding to S.hat.
#' @return A vector of expected values.
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))

  c(cbind(1, S.hat) %*% grid.diff)
}

# *** Estimators ***
estimate_grf = function(data, data.test) {
  num.trees = 500
  start = Sys.time()
  fit = survival_forest(data$X, data$Y, data$D, num.trees = num.trees)
  pp = predict(fit, data.test$X)
  pp.na = predict(fit, data.test$X, prediction.type = "Nelson-Aalen")
  end = Sys.time()

  elapsed.sec = difftime(end, start, units = "secs")
  ehat = expected_survival(pp$predictions, pp$failure.times)
  df.conc = data.frame(Y = data.test$Y, D = data.test$D, M = rowSums(-log(pp.na$predictions)))
  err = 1 - concordance(Surv(Y, D) ~ M, data = df.conc, reverse = TRUE)$concordance

  list(ehat = ehat, elapsed.sec = elapsed.sec, err = err)
}

estimate_SRC = function(data, data.test) {
  ntree = 500
  df = data.frame(Y=data$Y, D=data$D, x=data$X)
  df.test = data.frame(x=data.test$X)
  start = Sys.time()
  fit = rfsrc(Surv(Y, D) ~ ., ntree = ntree, nsplit = 0, data = df)
  pp = predict(fit, df.test)
  end = Sys.time()

  elapsed.sec = difftime(end, start, units = "secs")
  ehat = expected_survival(pp$survival, pp$time.interest)
  df.conc = data.frame(Y = data.test$Y, D = data.test$D, M = pp$predicted)
  err = 1 - concordance(Surv(Y, D) ~ M, data = df.conc, reverse = TRUE)$concordance

  list(ehat = ehat, elapsed.sec = elapsed.sec, err = err)
}

estimate_ranger = function(data, data.test) {
  num.trees = 500
  df = data.frame(Y=data$Y, D=data$D, x=data$X)
  df.test = data.frame(x=data.test$X)
  start = Sys.time()
  fit = ranger(Surv(Y, D) ~ ., data = df, num.trees = num.trees)
  pp = predict(fit, df.test)
  end = Sys.time()

  elapsed.sec = difftime(end, start, units = "secs")
  ehat = expected_survival(pp$survival, timepoints(fit))
  df.conc = data.frame(Y = data.test$Y, D = data.test$D, M = rowSums(pp$chf))
  err = 1 - concordance(Surv(Y, D) ~ M, data = df.conc, reverse = TRUE)$concordance

  list(ehat = ehat, elapsed.sec = elapsed.sec, err = err)
}

estimators = list(grf = estimate_grf,
                  SRC = estimate_SRC,
                  ranger = estimate_ranger)
# *** bench ***
grid = expand.grid(
  n = c(2000),
  p = c(5),
  n.test = c(2000),
  dgp = c("type1", "type2", "type3", "type4"),
  estimator = names(estimators),
  stringsAsFactors = FALSE
)
print(grid)

out = list()
n.sim = 150
print(Sys.time())
for (i in 1:nrow(grid)) {
  print(paste("grid:", i))
  n = grid$n[i]
  p = grid$p[i]
  n.test = grid$n.test[i]
  dgp = grid$dgp[i]
  estimator = grid$estimator[i]

  for (sim in 1:n.sim) {
    print(paste("sim", sim))
    data = generate_survival_data(n, p, dgp = dgp, n.mc = 10)
    data.test = generate_survival_data(n.test, p, dgp = dgp, n.mc = 1e5)
    est = estimators[[estimator]](data, data.test)

    df = data.frame(mse = mean((est$ehat - data.test$ET)^2),
                    err = est$err,
                    elapsed.sec = as.numeric(est$elapsed.sec),
                    n = n,
                    p = p,
                    n.test = n.test,
                    dgp = dgp,
                    estimator = estimator,
                    sim = sim,
                    stringsAsFactors = FALSE)
    out = c(out, list(df))
  }
}
print(Sys.time())
out.df = do.call(rbind, out)

write.csv(out.df, gzfile("bench_survival.csv.gz"), row.names = FALSE)
