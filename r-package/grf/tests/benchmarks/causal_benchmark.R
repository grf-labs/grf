# grf benchmark script
# This script benchmarks timings for:
# - regression forest
# - causal forest
# - survival forest
# - quantile forest
# - local linear regression forest
# This script also benchmarks statistical performance for:
# - regression forest
# - causal forest
# Usage: `(benchmarks)$ Rscript causal_benchmark.R`

library(grf)
set.seed(1)

reps <- 10
num.trees <- 2000

results.raw <- lapply(1:reps, function(iter) {
  data <- generate_causal_data(4000, 10, dgp = "simple")
  data.test <- generate_causal_data(4000, 10, dgp = "simple")

  # regression forest
  Y.time <- system.time(forest.Y <- regression_forest(data$X, data$Y,
                                                     num.trees = num.trees, num.threads=1))
  Y.time.pred <- system.time(Y.test <- predict(forest.Y, newdata = data.test$X, num.threads=1))
  Y.time.ci <- system.time(Y.test.ci <- predict(forest.Y, data.test$X,
                                               estimate.variance = TRUE, num.threads=1))
  Y.hat <- predict(forest.Y)$predictions

  Y.error <- c(
    oob = mean((Y.hat - data$m)^2),
    test = mean((Y.test$predictions - data.test$m)^2),
    stdratio = mean((Y.test.ci$predictions - data.test$m)^2 /
                      Y.test.ci$variance.estimates)
  )

  # causal forest
  tau.time <- system.time(forest.tau <- causal_forest(data$X, data$Y, data$W,
                                                     Y.hat=Y.hat, W.hat = data$e,
                                                     num.trees = num.trees, num.threads=1))
  tau.time.pred <- system.time(tau.test <- predict(forest.tau, data.test$X, num.threads=1))
  tau.time.ci <- system.time(tau.test.ci <- predict(forest.tau, data.test$X,
                                                   estimate.variance = TRUE, num.threads=1))

  tau.error <- c(
    oob = mean((predict(forest.tau)$predictions - data$tau)^2),
    test = mean((tau.test$predictions - data.test$tau)^2),
    stdratio = mean((tau.test.ci$predictions - data.test$tau)^2 /
                      tau.test.ci$variance.estimates)
  )

  # survival forest
  survival.time <- system.time(sf <- survival_forest(data$X[1:1500, ], data$Y[1:1500], data$W[1:1500],
                                                     num.trees = num.trees, num.threads=1))
  survival.time.pred <- system.time(sf.test <- predict(sf, newdata = data.test$X, num.threads=1))

  # quantile forest
  quantile.time <- system.time(qf <- quantile_forest(data$X[1:3000, ], data$Y[1:3000],
                                                     num.trees = num.trees, num.threads=1))
  quantile.time.pred <- system.time(qf.test <- predict(qf, newdata = data.test$X, num.threads=1))

  # local linear regression forest
  ll.Y.time <- system.time(ll.forest.Y <- ll_regression_forest(data$X, data$Y,
                                                               num.trees = num.trees, num.threads=1))
  ll.Y.time.pred <- system.time(ll.Y.test <- predict(ll.forest.Y, newdata = data.test$X[1:1000,], num.threads=1))

  error <- rbind(Y.error, tau.error)
  time <- rbind(Y.time, Y.time.pred, Y.time.ci,
                tau.time, tau.time.pred, tau.time.ci,
                survival.time, survival.time.pred,
                quantile.time, quantile.time.pred,
                ll.Y.time, ll.Y.time.pred)

  print("done with a rep!")
  list(error, time[,1:3])
})

avg.error = Reduce(function(a, b) a + b, lapply(results.raw, function(rr) rr[[1]])) / reps
avg.time = Reduce(function(a, b) a + b, lapply(results.raw, function(rr) rr[[2]])) / reps

se.error = sqrt((Reduce(function(a, b) a + b,
                        lapply(results.raw, function(rr) (rr[[1]] - avg.error)^2)) /
                   (reps * (reps - 1))))
se.time = sqrt((Reduce(function(a, b) a + b,
                       lapply(results.raw, function(rr) (rr[[2]] - avg.time)^2)) /
                  (reps * (reps - 1))))

err.print = outer(1:nrow(avg.error), 1:ncol(avg.error), Vectorize(function(ii, jj) {
  paste(sprintf("%.3f", avg.error[ii, jj]), "+/-", sprintf("%.3f", se.error[ii, jj]))
}))

time.print = outer(1:nrow(avg.time), 1:ncol(avg.time), Vectorize(function(ii, jj) {
  paste(sprintf("%.3f", avg.time[ii, jj]), "+/-", sprintf("%.3f", se.time[ii, jj]))
}))

colnames(err.print) = colnames(avg.error)
rownames(err.print) = rownames(avg.error)
colnames(time.print) = colnames(avg.time)
rownames(time.print) = rownames(avg.time)

err.print
time.print
