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
D <- as.integer(failure.time <= censor.time)
failure.times <- sort(unique(Y[D == 1]))
Y.relabeled <- findInterval(Y, failure.times)

sfit <- survfit(Surv(Y, D) ~ 1, data = data.frame(Y, D))

length(failure.times) # num_failures
dput(c(Y.relabeled, D)) # data_matrix
dput(summary(sfit)$surv) # expected_predictions Kaplan-Meier
dput(exp(-sfit$cumhaz[D[order(Y)] == 1])) # expected_predictions Nelson-Aalen
```
