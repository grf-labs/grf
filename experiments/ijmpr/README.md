_This folder has replication files for the paper "Estimating Treatment Effect Heterogeneity using Causal Forests: An Application to Stress Resilience" by Sverdrup, Petukhova, and Wager._

The file `analysis.R` contains example code to run the type of analysis described in the paper using synthetic data. This script relies on the packages `"grf", "maq", "ggplot2"`.

The file `synthetic_data.csv` was generated using the following code and is not intended to bear resemblance to the Army STARRS-LS data used in the paper.

```R
n = 4000
X = cbind(round(matrix(rnorm(n * 5), n, 5), 2), matrix(rbinom(n * 4, 1, 0.5), n, 4))
W = rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
Y = 1 - rbinom(n, 1, 1 / (1 + exp((pmax(2 * X[, 1], 0) * W + 1))))
colnames(X) = make.names(1:ncol(X))
write.csv(cbind(outcome = Y, treatment = W, X), "synthetic_data.csv", row.names = FALSE)
```
