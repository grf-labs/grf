# Datasets

The datasets in the `resources` directory were generated through the following R scripts:


```
# quantile_data.csv
p = 10
n = 1000
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
Y = (2 * rbinom(n, 1, 0.5) - 1) * (1 + 9 * (X[,1] > 0))
```

```
# quantile_data_MIA.csv
p <- 10
n <- 1000
X <- matrix(round(2 * runif(n * p) - 1, 2), n, p)
Y <- (2 * rbinom(n, 1, 0.5) - 1) * (1 + 9 * (X[,1] > 0))
nmissing <- 500
X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN
```

```
# probability_data.csv
p <- 10
n <- 1000
X <- matrix(round(rnorm(n * p), 2), n, p)
Y <- rbinom(n, 5, exp(-abs(X[, 1])))
```

```
# causal_data.csv
p = 10
n = 1000
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
W = rbinom(n, 1, 1/(1 + exp(X[,3])))
Y = round(100 * (2*W - 1) * (X[,1] > 0) + X[,2] + rnorm(n), 2)
```

```
# causal_data_MIA.csv
p <- 10
n <- 1000
X <- matrix(round(2 * runif(n * p) - 1, 2), n, p)
W <- rbinom(n, 1, 1/(1 + exp(X[,3])))
Y <- round(100 * (2*W - 1) * (X[,1] > 0) + X[,2] + rnorm(n), 2)
nmissing <- 500
X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN
```

```
# multi_causal_data.csv
n <- 500
p <- 5
X <- matrix(round(rnorm(n * p), 1), n, p)
W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
Y <- round(X[, 1] + X[, 2] * (W == "B") + X[, 3] * (W == "C") + runif(n), 1)
sample.weights <- sample(c(1, 5), n, TRUE)
W.matrix <- stats::model.matrix(~ W - 1)
data <- cbind(X, Y, W.matrix[, -1], sample.weights)
```

```
# gaussian_data.csv
n = 500
p = 11
X = matrix(rnorm(n * p), n , p)
```

```
# gaussian_data_shrunk.csv
n = 50
p = 11
X = matrix(rnorm(n * p), n , p)
```

```
# regression_data_MIA.csv
n <- 1000
p <- 5
X <- matrix(round(rnorm(n * p), 2), n, p)
Y <- round(X[, 1] * rnorm(n), 2)
nmissing <- 150
X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN
```

```
# survival_data.csv
# survival_data_MIA.csv
set.seed(123)
n <- 1000
p <- 5
X <- matrix(round(rnorm(n * p), 2), n, p)
failure.time <- -log(runif(n)) * exp(0.1 * X[, 1])
censor.time <- rexp(n)
Y <- round(pmin(failure.time, censor.time), 2)
D <- as.integer(failure.time <= censor.time)

nmissing <- 150
X[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN
```

```
# causal_survival_data.csv
# causal_survival_data_MIA.csv
set.seed(42)
n <- 500
p <- 5
dgp <- "simple1"
data <- generate_causal_survival_data(n = n, p = p, dgp = dgp)
X <- round(data$X, 2)
X.na <- X
Y <- round(data$Y, 2)
W <- data$W
D <- data$D
cs.forest <- causal_survival_forest(X, Y, W, D, num.trees = 50, seed = 42, num.threads = 4)
eta <- cs.forest$eta
nmissing <- 150
X.na[cbind(sample(1:n, nmissing), sample(1:p, nmissing, replace = TRUE))] <- NaN
cs.forest.na <- causal_survival_forest(X.na, Y, W, D, num.trees = 50, seed = 42, num.threads = 4)
eta.na <- cs.forest.na$eta

write.table(cbind(X, Y, W, D, eta[[1]], eta[[2]]), row.names = FALSE, col.names = FALSE,
            file = "causal_survival_data.csv")
write.table(cbind(X.na, Y, W, D, eta.na[[1]], eta.na[[2]]), row.names = FALSE, col.names = FALSE,
            na = "NaN",
            file = "causal_survival_data_MIA.csv")
```
