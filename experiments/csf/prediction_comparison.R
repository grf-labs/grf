# Produce a 1 by 3 plot of estimator predictions vs ground truth for a given DGP
# To reproduce top panel of Figure 1 in paper set `dgp = "type2"`
# To reproduce bottom panel of Figure 1 in paper set `dgp = "type3"`

rm(list = ls())
library(grf)
library(randomForestSRC)
library(survival)
library(ggplot2)
set.seed(123)

# *** Comparison methods ***
source("comparison_estimators.R")
estimators = list(estimate_rfsrc_X_W = estimate_rfsrc_X_W,
                  estimate_rfsrc_twin = estimate_rfsrc_twin,
                  estimate_grf = estimate_grf)
estimator.names = c("SRC", "VT", "CSF")

out = list()
n = 2000
p = 5
n.test = 2000
dgp = "type2"
# dgp = "type3"
n.mc = 100000

data = generate_causal_survival_data(n = n, p = p, dgp = dgp, n.mc = 10)
data.test = generate_causal_survival_data(n = n.test, p = p, dgp = dgp, n.mc = n.mc)
true.cate = data.test$cate
for (j in 1:length(estimators)) {
  estimator.name = estimator.names[j]
  predictions = estimators[[j]](data, data.test)
  out[[j]] = data.frame(
    predictions = predictions,
    estimator.name = estimator.name,
    true.cate = true.cate,
    mse = mean((predictions - true.cate)^2),
    cor = cor(predictions, true.cate)
    )
}
out = do.call(rbind, out)
out$label = paste0(out$estimator.name,
                   " (mse: ",
                   round(100 * out$mse, 2), # scale by 100
                   ", cor: ",
                   round(out$cor, 2)
                   ,")"
                   )
# Same order as rest of paper
out$label = factor(out$label, levels = unique(out$label)[c(2, 1, 3)])

ggplot(out, aes(y = predictions, x = true.cate)) +
  geom_point(size = 0.1) +
  geom_abline(intercept = 0, slope = 1, col = "red", lty = 3) +
  facet_wrap(. ~ label, ncol = 3) +
  theme_bw() +
  xlab("True effect") +
  # xlab("") +
  ylab("Estimated effect")

ggsave(paste0("prediction_comparsion_", dgp, ".pdf"), width = 6, height = 3)
