rm(list = ls())
set.seed(42)
library(grf)
library(microbenchmark)
library(xtable)

get_data = function(n, p, M = c("20", "130", "260", "500")) {
  M = match.arg(M)
  X = round(matrix(runif(n * p), n, p), 5)
  if (M == 20) {
    failure.time = rpois(n, 10 * X[, 1])
    censor.time = round(50 * rexp(n), 0)
  }
  if (M == 130) {
    failure.time = rpois(n, 100 * X[, 1])
    censor.time = round(500 * rexp(n), 0)
  }
  if (M == 260) {
    failure.time = rpois(n, 220 * X[, 1])
    censor.time = round(1000 * rexp(n), 0)
  }
  if (M == 500) {
    failure.time = rpois(n, 450 * X[, 1])
    censor.time = round(2000 * rexp(n), 0)
  }

  Y = pmin(failure.time, censor.time)
  D = as.integer(failure.time <= censor.time)

  list(X = X, Y = Y, D = D)
}

# n = 100000
# p = 50
# M = 500
# data = get_data(n, p, as.character(M))
# X = data$X; Y=data$Y; D=data$D
# mean(D)
# length(unique(Y[D==1]))
# length(unique(X[,1]))

fbench = function(data, fast.logrank) {
  survival_forest(data$X, data$Y, data$D, fast.logrank = fast.logrank,
                  compute.oob.predictions = FALSE,
                  num.trees = 1, sample.fraction = 1, mtry = ncol(data$X))
}

grid = expand.grid(
  n = c(20000, 50000, 250000),
  p = c(25, 50, 100),
  M = c(20, 130, 260, 500)
)

times = 10
out = list()
for (i in 1:nrow(grid)) {
  n = grid[i, ]$n
  p = grid[i, ]$p
  M = grid[i, ]$M
  data = get_data(n = n, p = p, M = as.character(M))

  bench.exact = microbenchmark(fbench(data, FALSE), times = times, unit = "seconds")
  time.exact = summary(bench.exact)$mean
  bench.approx = microbenchmark(fbench(data, TRUE), times = times, unit = "seconds")
  time.approx = summary(bench.approx)$mean

  diff = time.exact - time.approx
  ratio = time.exact / time.approx
  df = data.frame(
    metric = c("exact(s)", "approx(s)", "difference(s)", "speedup.factor"),
    value = c(time.exact, time.approx, diff, ratio),
    n = n,
    p = p,
    M = M)
  out = c(out, list(df))
}
out.df = do.call(rbind, out)
write.csv(out.df, "timing.csv", row.names = FALSE)

tab.df = reshape(
  out.df,
  idvar = c("n", "p", "M"),
  timevar = "metric",
  direction = "wide"
)


tab.df$n = as.character(tab.df$n)
tab.df$p = as.character(tab.df$p)
tab.df$M = as.character(tab.df$M)
print(xtable(tab.df[-6]), include.rownames = FALSE)


data = get_data(100000, 50, "260")
print(m1 <- microbenchmark(survival_forest(data$X, data$Y, data$D, fast.logrank = FALSE, num.trees = 500, compute.oob.predictions = FALSE), times = 2, unit = "seconds"))
print(m2 <- microbenchmark(survival_forest(data$X, data$Y, data$D, fast.logrank = TRUE, num.trees = 500, compute.oob.predictions = FALSE), times = 2, unit = "seconds"))
summary(m1)$mean / summary(m2)$mean
