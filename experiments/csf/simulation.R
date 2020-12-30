rm(list = ls())
library(grf)
library(randomForestSRC)
library(survival)
set.seed(123)

# *** Comparison methods ***
source("comparison_estimators.R")
estimators = list(estimate_rfsrc_X_W = estimate_rfsrc_X_W,
                  estimate_rfsrc_XW_W = estimate_rfsrc_XW_W,
                  estimate_rfsrc_twin = estimate_rfsrc_twin,
                  estimate_IPCW_grf = estimate_IPCW_grf,
                  estimate_grf = estimate_grf)

# *** Setup ***
out = list()
n.sim = 200
n.mc = 100000
grid = expand.grid(n = c(500, 1000, 2000, 5000),
                   p = 5,
                   n.test = 2000,
                   dgp = c("type1", "type2", "type3", "type4"),
                   stringsAsFactors = FALSE)

for (i in 1:nrow(grid)) {
  print(paste("grid", i))
  n = grid$n[i]
  p = grid$p[i]
  n.test = grid$n.test[i]
  dgp = grid$dgp[i]

  for (sim in 1:n.sim) {
    print(paste("sim", sim))
    data = generate_survival_data(n = n, p = p, dgp = dgp, n.mc = 10)
    data.test = generate_survival_data(n = n.test, p = p, dgp = dgp, n.mc = n.mc)
    true.cate = data.test$cate
    true.cate.sign = data.test$cate.sign
    estimator.output = list()
    for (j in 1:length(estimators)) {
      estimator.name = names(estimators)[j]
      predictions = estimators[[estimator.name]](data, data.test)
      correct.classification = sign(predictions) == true.cate.sign
      dfj = data.frame(
        estimator.name = estimator.name,
        mse = mean((predictions - true.cate)^2),
        classif.rate = mean(correct.classification, na.rm = TRUE) # NA: to ignore X1 < 0.3 in DGP 4.
        )
      estimator.output[[j]] = dfj
    }
    df = do.call(rbind, estimator.output)
    df$n = n
    df$p = p
    df$n.test = n.test
    df$dgp = dgp
    df$sim = sim

    out = c(out, list(df))
  }
}
out.df = do.call(rbind, out)

write.csv(out.df, gzfile("simulation.csv.gz"), row.names = FALSE)
