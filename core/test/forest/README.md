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
