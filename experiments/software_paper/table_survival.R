rm(list = ls())
library(xtable)

df = read.csv("bench_survival.csv.gz")
head(df)
# Set order
df$estimator = factor(df$estimator, levels = c("SRC", "ranger", "grf"))

tab.agg = aggregate(list(mse = 100*df$mse, time = df$elapsed.sec, err = 100*df$err),
                    by = list(n = df$n,
                              p = df$p,
                              n.test = df$n.test,
                              dgp = df$dgp,
                              estimator = df$estimator),
                    FUN = mean)
tab.agg.se = aggregate(list(mse.se = 100*df$mse, time.se = df$elapsed.sec, err.se = 100*df$err),
                              by = list(n = df$n,
                                        p = df$p,
                                        n.test = df$n.test,
                                        dgp = df$dgp,
                                        estimator = df$estimator),
                              FUN = function(x) sd(x) / sqrt(length(x)))
tab.agg$mse = paste0(format(round(tab.agg$mse, 1), nsmall=1), " (", format(round(tab.agg.se$mse.se, 2), nsmall=2), ")")
tab.agg$time = paste0(format(round(tab.agg$time, 1), nsmall=1), " (", format(round(tab.agg.se$time.se, 2), nsmall=2), ")")
tab.agg$err = paste0(format(round(tab.agg$err, 1), nsmall=1), " (", format(round(tab.agg.se$err.se, 2), nsmall=2), ")")

print(tab.agg, digits = 2)
apply(df[c("n", "p", "n.test", "dgp", "estimator")], 2, unique)

tab.out = rbind(
  cbind(tab.agg[c("dgp", "estimator")], metric = "RMST MSE", value = tab.agg[["mse"]]),
  cbind(tab.agg[c("dgp", "estimator")], metric = "error (%)", value = tab.agg[["err"]]),
  cbind(tab.agg[c("dgp", "estimator")], metric = "time (sec.)", value = tab.agg[["time"]])
)
head(tab.out)
wide = reshape(tab.out, timevar = c("estimator"), idvar = c("dgp", "metric"), direction = "wide")

# Table 1 - RMST simulation
print(
  xtable(wide[order(wide$dgp, wide$metric), ]),
  include.rownames = FALSE, digits = 2
)
