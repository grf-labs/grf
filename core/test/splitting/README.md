# Datasets

The datasets in the `resources` directory were generated through the following R scripts:


```R
# survival_data_logrank.csv
# survival_data_logrank_expected.csv

# This script loops through all splits points for all variables
# and uses the `survdiff` function from the `survival` package to compute the logrank
# statistic at the given split. The largest logrank statistic subject to a minimum
# of one failure on both sides of the split is saved in `survival_data_logrank_expected.csv`.
# p = 500 gives 500 logrank statistics computed on random partitions to check the SurvivalSplittingRule against.

library(survival)

n = 100
p = 500
X = matrix(round(rnorm(n * p), 4), n, p)
Y = round(rexp(n), 2)
C = rbinom(n, 1, 0.75)

df = data.frame(X, Y, C) # survival_data_logrank.csv
out = matrix(0, p) # survival_data_logrank_expected.csv

num.failures = length(Y[C==1])
for (var in 1:p) {
  sorted.samples = order(df[, var])
  split.values = sort(unique(df[, var]))
  num.splits = length(split.values) - 1
  logranks = matrix(0, num.splits, 2)
  for (i in 1:num.splits) {
    sval = split.values[i]
    df$split.group = df[, var] <= sval
    lr = survdiff(Surv(Y, C) ~ split.group, data = df)$chisq
    logranks[i, 1] = lr
    subset = df[df$split.group, ]
    logranks[i, 2] = sum(subset["C"] == 1)
  }
  count.failures.left = logranks[, 2]
  count.failures.right = num.failures - count.failures.left
  out[var, 1] = max(logranks[count.failures.left >= 1 & count.failures.right >= 1, 1])
}
```
