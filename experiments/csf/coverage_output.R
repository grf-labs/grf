# Run `coverage.R` to produce `coverage.csv.gz`.
# Table 3 is produced below.

rm(list = ls())
library(xtable)
df = read.csv("coverage.csv.gz")
apply(df[c("n", "num.trees", "dgp", "X.test")], 2, unique)

tab = aggregate(list(coverage = df$coverage),
                by = list(dgp = df$dgp,
                          Xi = df$X.test,
                          n.train = df$n,
                          num.trees = df$num.trees),
                FUN = mean)

# Table 3
options(digits = 2)
xtable(xtabs(coverage ~ dgp + Xi, tab))

