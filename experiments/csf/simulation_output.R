# Tabulate simulation results
# Run `simulation.R` to produce `simulation.csv.gz`.
# Table 1, Table 2, appendix Figure 5 and Figure 6 are produced below.

library(xtable)
rm(list = ls())

df = read.csv("simulation.csv.gz", stringsAsFactors = FALSE)
apply(df[c("n", "p", "n.test", "dgp", "estimator.name")], 2, unique)

old.names = c("estimate_rfsrc_X_W",
              "estimate_rfsrc_XW_W",
              "estimate_rfsrc_twin",
              "estimate_IPCW_grf",
              "estimate_grf")
new.names = c("SRC1",
              "SRC2",
              "VT",
              "IPCW",
              "CSF")
df$estimator = new.names[match(df$estimator.name, old.names)]
df$estimator.name = NULL
df$classif.error = 1 - df$classif.rate
df$mse = df$mse * 100

df.min = aggregate(list(mse.min = df$mse),
                   by = list(n = df$n,
                             p = df$p,
                             n.test = df$n.test,
                             dgp = df$dgp,
                             sim = df$sim),
                    FUN = min)
DF = merge(df, df.min)
DF$mse.excess = DF$mse / DF$mse.min

# Reorder estimator names
DF$estimator = factor(DF$estimator, levels = c("VT", "SRC1", "SRC2", "IPCW", "CSF"))

# Tables
DF.table = aggregate(list(MSE = DF$mse,
                          MSE.excess = DF$mse.excess,
                          CLF = DF$classif.error),
                          by = list(estimator = DF$estimator,
                                      n = DF$n,
                                      p = DF$p,
                                      n.test = DF$n.test,
                                      dgp = DF$dgp),
                          FUN = mean)
# Look at n = 2000
DF.table = DF.table[DF.table$n == 2000, ]
tab.mse = DF.table[c("dgp", "estimator", "MSE", "MSE.excess")]
tab.mse = rbind(cbind(tab.mse[, 1:2], metric = "mse", value = tab.mse[, 3]),
                cbind(tab.mse[, 1:2], metric = "mse.excess", value = tab.mse[, 4])
)
tab.mse = reshape(tab.mse, timevar = c("estimator"), idvar = c("dgp", "metric"), direction = "wide")
tab.mse = tab.mse[order(tab.mse$dgp, tab.mse$metric), ]
# Table 1 - MSE
print(xtable(tab.mse), include.rownames = FALSE)

tab.clf = DF.table[c("dgp", "estimator", "CLF")]
tab.clf = xtabs(CLF ~ dgp + estimator, tab.clf)
# Table 2 - classification error
print(xtable(tab.clf))

# *** Appendix plots ***
library(ggplot2)
# Figure 5 - Classification error for each n
ggplot(DF, aes(x = estimator, y = classif.error)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Classification error") +
  xlab("Method")
ggsave("simulation_comparison_classif_allN.pdf", width = 9, height = 9)

# Figure 6 - MSE for each n
ggplot(DF, aes(x = estimator, y = mse.excess)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Excess MSE") +
  xlab("Method")
ggsave("simulation_comparison_mse_allN.pdf", width = 9, height = 9)
