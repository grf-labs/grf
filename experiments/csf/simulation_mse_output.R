# Tabulate simulation results
# Run `simulation_mse.R` to produce `simulation.csv.gz`.
# Table 1, Table 2, appendix Figure 5 and Figure 6 are produced below.

library(xtable)
rm(list = ls())

df = read.csv("simulation.csv.gz", stringsAsFactors = FALSE)
apply(df[c("n", "p", "n.test", "dgp", "estimator", "target")], 2, unique)

df$CLF = 1 - df$classif.rate # classification error
df$MSE = df$MSE * 100

df.min = aggregate(list(MSE.min = df$MSE),
                   by = list(target = df$target,
                             n = df$n,
                             p = df$p,
                             n.test = df$n.test,
                             dgp = df$dgp,
                             sim = df$sim),
                    FUN = min)
DF = merge(df, df.min)
DF$MSE.excess = DF$MSE / DF$MSE.min

# Reorder estimator names
DF$estimator = factor(DF$estimator, levels = c("VT", "SRC1", "SRC2", "IPCW", "CSF"))

# ***  Tables ***
DF.table = aggregate(list(MSE = DF$MSE,
                          MSE.excess = DF$MSE.excess,
                          CLF = DF$CLF),
                     by = list(estimator = DF$estimator,
                               target = DF$target,
                               n = DF$n,
                               p = DF$p,
                               n.test = DF$n.test,
                               dgp = DF$dgp),
                     FUN = mean)

metrics.col = which(names(DF.table) %in% c("MSE", "MSE.excess", "CLF"))
DF.table.long = reshape(DF.table,
                        direction = "long",
                        varying = list(names(DF.table)[metrics.col]),
                        v.names = "val",
                        idvar = names(DF.table)[-metrics.col],
                        times = names(DF.table)[metrics.col])
rownames(DF.table.long) = NULL
colnames(DF.table.long)[colnames(DF.table.long) == "time"] = "metric"

# Look at n = 2000
# DF.table.long = subset(DF.table.long, n == 2000)
xtabs(val ~ dgp + estimator + metric + target, DF.table.long, subset = n == 2000)

# tab.mse = DF.table[c("dgp", "estimator", "MSE", "MSE.excess")]
# tab.mse = rbind(cbind(tab.mse[, 1:2], metric = "mse", value = tab.mse[, 3]),
                # cbind(tab.mse[, 1:2], metric = "mse.excess", value = tab.mse[, 4])
# )
# tab.mse = reshape(tab.mse, timevar = c("estimator"), idvar = c("dgp", "metric"), direction = "wide")
# tab.mse = tab.mse[order(tab.mse$dgp, tab.mse$metric), ]
# Table 1 - MSE
# print(xtable(tab.mse), include.rownames = FALSE)

# tab.clf = DF.table[c("dgp", "estimator", "CLF")]
# tab.clf = xtabs(CLF ~ dgp + estimator, tab.clf)
# Table 2 - classification error
# print(xtable(tab.clf))

# *** Appendix plots ***
library(ggplot2)
# Figure 5 - Classification error for each n
ggplot(subset(DF, target == "RMST"), aes(x = estimator, y = CLF)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Classification error") +
  xlab("Method")
ggsave("simulation_comparison_CLF.pdf", width = 9, height = 9)

# Figure 6 - excess MSE for each n
ggplot(subset(DF, target == "RMST"), aes(x = estimator, y = MSE.excess)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Excess MSE") +
  xlab("Method")
ggsave("simulation_comparison_MSEexcess.pdf", width = 9, height = 9)
