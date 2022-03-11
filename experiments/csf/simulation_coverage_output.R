# Run `simulation_coverage.R` to produce `coverage.csv.gz`.
# 95 % CI coverage table is produced below.

rm(list = ls())
library(xtable)
df = read.csv("coverage.csv.gz")
apply(df[c("target", "n", "p", "rho", "num.trees", "dgp", "X.test")], 2, unique)

tab = aggregate(list(coverage = df$coverage, width = df$width),
                by = list(target = df$target,
                          dgp = df$dgp,
                          Xi = df$X.test,
                          p = df$p,
                          rho = df$rho,
                          n.train = df$n,
                          num.trees = df$num.trees),
                FUN = mean)

# Table coverage and CI length
# RMST
print(xtable(
  cbind(xtabs(coverage ~ dgp + Xi, tab, subset = target == "RMST" & rho == 0),
      xtabs(width ~ dgp + Xi, tab, subset = target == "RMST" & rho == 0))
))
# SP
print(xtable(
  cbind(xtabs(coverage ~ dgp + Xi, tab, subset = target == "survival.probability" & rho == 0),
        xtabs(width ~ dgp + Xi, tab, subset = target == "survival.probability" & rho == 0))
))

# w correlated X's
# RMST
print(xtable(
  cbind(xtabs(coverage ~ dgp + Xi, tab, subset = target == "RMST" & rho == 0.5),
        xtabs(width ~ dgp + Xi, tab, subset = target == "RMST" & rho == 0.5))
))
# SP
print(xtable(
  cbind(xtabs(coverage ~ dgp + Xi, tab, subset = target == "survival.probability" & rho == 0.5),
        xtabs(width ~ dgp + Xi, tab, subset = target == "survival.probability" & rho == 0.5))
))
