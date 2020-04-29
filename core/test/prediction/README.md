# Datasets

The datasets in the test cases were generated through the following R scripts:


```R
# SurvivalPredictionStrategyTest.cpp

library(survival)
n <- 50
p <- 1
X <- matrix(rnorm(n * p), n, p)
failure.time <- -log(runif(n)) * exp(0.1 * X[, 1])
censor.time <- rexp(n)
Y <- pmin(failure.time, censor.time)
C <- as.integer(failure.time <= censor.time)
failure.times <- sort(unique(Y[C == 1]))
Y.relabeled <- findInterval(Y, failure.times, rightmost.closed = FALSE, all.inside = FALSE, left.open = FALSE)

kaplan.meier <- survfit(Surv(Y, C) ~ 1, data = data.frame(Y, C))

length(failure.times) # num_failures
dput(c(Y.relabeled, C)) # data_matrix
dput(summary(kaplan.meier)$surv) # expected_predictions
```
