# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
rm(list = ls())
set.seed(42)
library(randomForestSRC) # 2.9.3
library(ranger) # 0.12.1
# library(tuneRanger) # 0.5 # Ignore, does not install w/o error on a linux box
library(grf) # 1.2.0.0

# *** Estimators ***
estimate_grf = function(data, data.test, honesty = TRUE, tune = FALSE) {
  num.trees = 500
  tune.parameters = if (tune) "all" else "none"
  start = Sys.time()
  fit = regression_forest(data$X, data$Y, num.trees = num.trees, honesty = honesty, tune.parameters = tune.parameters)
  pp = predict(fit, data.test$X)
  end = Sys.time()
  elapsed.sec = difftime(end, start, units = "secs")

  list(y.hat = pp$predictions, elapsed.sec = elapsed.sec)
}

estimate_SRC = function(data, data.test) {
  ntree = 500
  df = data.frame(Y=data$Y, x=data$X)
  df.test = data.frame(x=data.test$X)
  start = Sys.time()
  fit = rfsrc(Y ~ ., ntree = ntree, nsplit = 0, data = df)
  pp = predict(fit, df.test)
  end = Sys.time()
  elapsed.sec = difftime(end, start, units = "secs")

  list(y.hat = pp$predicted, elapsed.sec = elapsed.sec)
}

estimate_ranger = function(data, data.test, tune = FALSE) {
  num.trees = 500
  df = data.frame(Y=data$Y, x=data$X)
  df.test = data.frame(x=data.test$X)
  start = Sys.time()
  if (tune) {
    task = makeRegrTask(data = df, target = "Y")
    res <- tryCatch({
      tuneRanger(task, build.final.model = FALSE)
      },
      error = function(e) {
        warning(paste0("ranger tuning threw the following error (fitting with default values): \n", e))
        return(NULL)
    })
    if (is.null(res)) {
      fit = ranger(Y ~ ., data = df, num.trees = num.trees)
      } else {
      fit = ranger(Y ~ ., data = df, num.trees = num.trees, mtry = res$recommended.pars$mtry,
                   sample.fraction = res$recommended.pars$sample.fraction)
    }
  } else {
    fit = ranger(Y ~ ., data = df, num.trees = num.trees)
  }
  pp = predict(fit, df.test)
  end = Sys.time()
  elapsed.sec = difftime(end, start, units = "secs")

  list(y.hat = pp$predictions, elapsed.sec = elapsed.sec)
}

estimators = list(grf = estimate_grf,
                  grf.adaptive = function (data, data.test) estimate_grf(data, data.test, honesty = FALSE),
                  grf.tuned = function (data, data.test) estimate_grf(data, data.test, tune = TRUE),
                  SRC = estimate_SRC,
                  ranger = estimate_ranger
                  )
# *** bench ***
grid = expand.grid(
  n = c(10000),
  p = c(30),
  n.test = c(2000),
  dgp = c("kunzel", "ai1", "ai2"),
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
    data = generate_causal_data(n, p, dgp = dgp)
    data.test = generate_causal_data(n.test, p, dgp = dgp)
    est = estimators[[estimator]](data, data.test)

    df = data.frame(mse = mean((est$y.hat - data.test$Y)^2),
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

write.csv(out.df, gzfile("bench_regression.csv.gz"), row.names = FALSE)
