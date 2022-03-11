# Tabulate simulation results
# Run `simulation_mse.R` to produce `simulation.csv.gz`.
# MSE/classification error tables + boxplot figures are produced below.

library(xtable)
rm(list = ls())

df = read.csv("simulation.csv.gz", stringsAsFactors = FALSE)
apply(df[c("n", "p", "rho", "n.test", "dgp", "estimator", "target")], 2, unique)

df$CLF = 1 - df$classif.rate # classification error
df$MSE = df$MSE * 100

df.min = aggregate(list(MSE.min = df$MSE),
                   by = list(target = df$target,
                             n = df$n,
                             p = df$p,
                             rho = df$rho,
                             n.test = df$n.test,
                             dgp = df$dgp,
                             sim = df$sim),
                    FUN = min)
DF = merge(df, df.min)
DF$MSE.excess = DF$MSE / DF$MSE.min

# Reorder estimator names
DF$estimator = factor(DF$estimator, levels = c("VT", "SRC1", "SRC2", "IPCW", "CSF"))

DF.table = aggregate(list(MSE = DF$MSE,
                          MSE.excess = DF$MSE.excess,
                          CLF = DF$CLF),
                     by = list(estimator = DF$estimator,
                               target = DF$target,
                               n = DF$n,
                               p = DF$p,
                               rho = DF$rho,
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

# ***  Tables ***
DF.table.long = subset(DF.table.long, n == 2000 & rho == 0)

# Table MSE + Excess MSE
t1 = subset(DF.table.long, metric != "CLF" & target == "RMST")
t1w = reshape(t1[c("estimator", "dgp", "metric", "val")],
              timevar = "estimator", idvar = c("dgp", "metric"), direction = "wide")
print(xtable(
  t1w[order(t1w$dgp, t1w$metric), ]),
  include.rownames = FALSE)

t2 = subset(DF.table.long, metric != "CLF" & target == "survival.probability")
t2w = reshape(t2[c("estimator", "dgp", "metric", "val")],
              timevar = "estimator", idvar = c("dgp", "metric"), direction = "wide")
print(xtable(
  t2w[order(t2w$dgp, t2w$metric), ]),
  include.rownames = FALSE)

# Table classification error
print(xtable(
    xtabs(val ~ dgp + estimator,
          DF.table.long,
          subset = metric == "CLF" & target == "RMST")))
print(xtable(
    xtabs(val ~ dgp + estimator,
          DF.table.long,
          subset = metric == "CLF" & target == "survival.probability")))

# *** Appendix plots ***
library(ggplot2)
## Independent covariates
ggplot(subset(DF, target == "RMST" & rho == 0), aes(x = estimator, y = CLF)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Classification error (RMST)") +
  xlab("Method")
ggsave("sim_CLF_RMST.pdf", width = 9, height = 9)

ggplot(subset(DF, target == "RMST" & rho == 0), aes(x = estimator, y = MSE.excess)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Excess MSE (RMST)") +
  xlab("Method")
ggsave("sim_MSEexcess_RMST.pdf", width = 9, height = 9)

ggplot(subset(DF, target == "survival.probability" & rho == 0), aes(x = estimator, y = CLF)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Classification error (survival probabilities)") +
  xlab("Method")
ggsave("sim_CLF_SP.pdf", width = 9, height = 9)

ggplot(subset(DF, target == "survival.probability" & rho == 0), aes(x = estimator, y = MSE.excess)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Excess MSE (survival probabilities)") +
  xlab("Method")
ggsave("sim_MSEexcess_SP.pdf", width = 9, height = 9)

## Correlated covariates
ggplot(subset(DF, target == "RMST" & rho != 0), aes(x = estimator, y = CLF)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Classification error (RMST)") +
  xlab("Method")
ggsave("sim_CLF_RMST_cor.pdf", width = 9, height = 9)

ggplot(subset(DF, target == "RMST" & rho != 0), aes(x = estimator, y = MSE.excess)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Excess MSE (RMST)") +
  xlab("Method")
ggsave("sim_MSEexcess_RMST_cor.pdf", width = 9, height = 9)

ggplot(subset(DF, target == "survival.probability" & rho != 0), aes(x = estimator, y = CLF)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Classification error (survival probabilities)") +
  xlab("Method")
ggsave("sim_CLF_SP_cor.pdf", width = 9, height = 9)

ggplot(subset(DF, target == "survival.probability" & rho != 0), aes(x = estimator, y = MSE.excess)) +
  geom_boxplot() +
  facet_wrap(dgp ~ n, scales = "free_y") +
  theme_bw() +
  ylab("Excess MSE (survival probabilities)") +
  xlab("Method")
ggsave("sim_MSEexcess_SP_cor.pdf", width = 9, height = 9)
