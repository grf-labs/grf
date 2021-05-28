rm(list = ls())
library(xtable)

df = read.csv("bench_regression.csv.gz")
head(df)
# Set order
df$estimator = factor(df$estimator, levels = c("SRC", "ranger", "grf", "grf.adaptive", "grf.tuned"))

tab.agg = aggregate(list(mse = df$mse, time = df$elapsed.sec),
                    by = list(n = df$n,
                              p = df$p,
                              n.test = df$n.test,
                              dgp = df$dgp,
                              estimator = df$estimator),
                    FUN = mean)
tab.agg.se = aggregate(list(mse.se = df$mse, time.se = df$elapsed.sec),
                             by = list(n = df$n,
                                       p = df$p,
                                       n.test = df$n.test,
                                       dgp = df$dgp,
                                       estimator = df$estimator),
                             FUN = function(x) sd(x) / sqrt(length(x)))
tab.agg$mse = paste0(format(round(tab.agg$mse, 1), nsmall=1), " (", format(round(tab.agg.se$mse.se, 2), nsmall=2), ")")
tab.agg$time = paste0(format(round(tab.agg$time, 1), nsmall=1), " (", format(round(tab.agg.se$time.se, 2), nsmall=2), ")")

print(tab.agg, digits = 2)
apply(df[c("n", "p", "n.test", "dgp", "estimator")], 2, unique)

tab.out = rbind(
  cbind(tab.agg[c("dgp", "estimator")], metric = "MSE", value = tab.agg[["mse"]]),
  cbind(tab.agg[c("dgp", "estimator")], metric = "time (sec.)", value = tab.agg[["time"]])
)
head(tab.out)
wide = reshape(tab.out, timevar = c("estimator"), idvar = c("dgp", "metric"), direction = "wide")

# Table 2 - MSE simulation
print(
  xtable(wide[order(wide$dgp, wide$metric), ]),
  include.rownames = FALSE, digits = 2
)
