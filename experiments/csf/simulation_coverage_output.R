# Run `simulation_coverage.R` to produce `coverage.csv.gz`.
# Table 3 is produced below.

rm(list = ls())
library(xtable)
df = read.csv("coverage.csv.gz")
apply(df[c("target", "n", "num.trees", "dgp", "X.test")], 2, unique)

tab = aggregate(list(coverage = df$coverage),
                by = list(target = df$target,
                          dgp = df$dgp,
                          Xi = df$X.test,
                          n.train = df$n,
                          num.trees = df$num.trees),
                FUN = mean)

# Table 3
options(digits = 2)
xtabs(coverage ~ dgp + Xi + target, tab)
