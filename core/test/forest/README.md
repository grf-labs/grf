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
# causal_data.csv
p = 10
n = 1000
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
W = rbinom(n, 1, 1/(1 + exp(X[,3])))
Y = round(100 * (2*W - 1) * (X[,1] > 0) + X[,2] + rnorm(n), 2)
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
