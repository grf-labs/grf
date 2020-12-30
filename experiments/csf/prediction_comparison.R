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

data = generate_survival_data(n = n, p = p, dgp = dgp, n.mc = 10)
data.test = generate_survival_data(n = n.test, p = p, dgp = dgp, n.mc = n.mc)
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

# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
#
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
#
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] ggplot2_3.3.2         survival_3.2-3        randomForestSRC_2.9.3 grf_1.2.0.0
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.5       rstudioapi_0.11  magrittr_1.5     splines_3.6.1    tidyselect_1.1.0 munsell_0.5.0
# [7] colorspace_1.4-1 lattice_0.20-38  R6_2.4.0         rlang_0.4.7      dplyr_0.8.3      tools_3.6.1
# [13] parallel_3.6.1   grid_3.6.1       gtable_0.3.0     withr_2.1.2      digest_0.6.20    assertthat_0.2.1
# [19] tibble_2.1.3     crayon_1.3.4     Matrix_1.2-17    purrr_0.3.3      glue_1.3.1       labeling_0.3
# [25] compiler_3.6.1   pillar_1.4.2     scales_1.0.0     pkgconfig_2.0.2
