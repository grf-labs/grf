# Datasets

The datasets is this directory were generated through the following R scripts:

```
# causal_survival_oob_predictions.csv
# causal_survival_predictions.csv
set.seed(42)
n <- 500
p <- 5
dgp <- "simple1"
data <- generate_survival_data(n = n, p = p, dgp = dgp)
cs.forest <- causal_survival_forest(round(data$X, 2), round(data$Y, 2), data$W, data$D,
                                    num.trees = 50, seed = 42, num.threads = 4)


write.table(predict(cs.forest)$predictions, file = "causal_survival_oob_predictions.csv", row.names = FALSE, col.names = FALSE)
write.table(predict(cs.forest, round(data$X, 2))$predictions, file = "causal_survival_predictions.csv", row.names = FALSE, col.names = FALSE)
```

```
# causal_survival_oob_predictions_grid.csv
# causal_survival_predictions_grid.csv
set.seed(42)
n <- 500
p <- 5
dgp <- "simple1"
data <- generate_survival_data(n = n, p = p, dgp = dgp)
failure.times <- seq(min(data$Y), max(data$Y), length.out = 5)
cs.forest.grid <- causal_survival_forest(round(data$X, 2), round(data$Y, 2), data$W, data$D,
                                         failure.times = failure.times,
                                         num.trees = 50, seed = 42, num.threads = 4)

write.table(predict(cs.forest.grid)$predictions, file = "causal_survival_oob_predictions_grid.csv", row.names = FALSE, col.names = FALSE)
write.table(predict(cs.forest.grid, round(data$X, 2))$predictions, file = "causal_survival_predictions_grid.csv", row.names = FALSE, col.names = FALSE)
```
