# Produce a 1 by 3 plot of estimator predictions vs ground truth for a given DGP
# To reproduce top panel of Figure 1 in paper set `dgp = "type2"`
# To reproduce bottom panel of Figure 1 in paper set `dgp = "type3"`

rm(list = ls())
library(grf)
library(randomForestSRC)
library(ggplot2)
set.seed(123)

# *** Comparison methods ***
source("estimators.R")
estimators = list(SRC = SRC1,
                  VT = VT,
                  CSF = CSF)

out = list()
n = 5000
p = 15
n.test = 2000
dgp = "type2"
# dgp = "type3"

data = generate_causal_survival_data(n = n, p = p, dgp = dgp, n.mc = 1)
data$Y = round(data$Y, 2)
data.test = generate_causal_survival_data(n = n.test, p = p, dgp = dgp, n.mc = 100000)
true.cate = data.test$cate
for (j in 1:length(estimators)) {
  estimator = names(estimators)[j]
  predictions = estimators[[j]](data, data.test)$pp
  out[[j]] = data.frame(
    predictions = predictions,
    estimator = estimator,
    true.cate = true.cate,
    mse = mean((predictions - true.cate)^2),
    cor = cor(predictions, true.cate)
    )
}
out = do.call(rbind, out)
out$label = paste0(out$estimator,
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
  geom_abline(intercept = 0, slope = 1, col = "red", lty = 1) +
  facet_wrap(. ~ label, ncol = 3) +
  theme_bw() +
  # xlab("True effect") +
  xlab("") +
  ylab("Estimated effect")

ggsave(paste0("prediction_comparsion_", dgp, ".pdf"), width = 6, height = 3)
