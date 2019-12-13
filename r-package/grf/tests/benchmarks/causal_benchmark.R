# grf benchmark script
# Usage: `(benchmarks)$ Rscript causal_benchmark.R`

library(grf)
source("../../../../experiments/benchmarking/dgps.R")
set.seed(1)

reps <- 10
num.trees <- 2000

results.raw <- lapply(1:reps, function(iter) {
  data <- gen_data(4000, 10, dgp = "simple")
  data.test <- gen_data(4000, 10, dgp = "simple")

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

  error <- rbind(Y.error, tau.error)
  time <- rbind(Y.time, Y.time.pred, Y.time.ci,
                tau.time, tau.time.pred, tau.time.ci)

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
