rm(list = ls())
library(grf)
library(randomForestSRC)
set.seed(123)

# *** Comparison methods ***
source("estimators.R")
estimators = list(SRC1 = SRC1,
                  SRC2 = SRC2,
                  VT = VT,
                  IPCW = IPCW,
                  CSF = CSF)

# *** Setup ***
out = list()
n.sim = 250
grid = expand.grid(n = c(500, 1000, 2000, 5000),
                   p = 15,
                   rho = c(0, 0.5),
                   n.test = 2000,
                   dgp = c("type1", "type2", "type3", "type4"),
                   stringsAsFactors = FALSE)

for (i in 1:nrow(grid)) {
  print(Sys.time())
  print(paste("grid", i, "of", nrow(grid)))
  print(grid[i, ])
  n = grid$n[i]
  p = grid$p[i]
  n.test = grid$n.test[i]
  dgp = grid$dgp[i]
  rho = grid$rho[i]

  data.test = generate_causal_survival_data(n = n.test, p = p, dgp = dgp, rho = rho, n.mc = 100000)
  true.cate = data.test$cate
  true.cate.prob = data.test$cate.prob
  true.cate.sign = data.test$cate.sign
  for (sim in 1:n.sim) {
    print(paste("sim", sim))
    data = generate_causal_survival_data(n = n, p = p, dgp = dgp, rho = rho, n.mc = 1)
    data$Y = round(data$Y, 2)
    estimator.output = list()
    for (j in 1:length(estimators)) {
      estimator = names(estimators)[j]
      predictions = estimators[[j]](data, data.test)
      dfj = data.frame(
        estimator = estimator,
        target = c("RMST", "survival.probability"),
        MSE = c(mean((predictions$pp - true.cate)^2),
                mean((predictions$pp.prob - true.cate.prob)^2)),
        classif.rate = c(mean(sign(predictions$pp) == true.cate.sign, na.rm = TRUE), # na.rm: to ignore X1 < 0.3 in DGP 4.
                         mean(sign(predictions$pp.prob) == true.cate.sign, na.rm = TRUE))
        )
      estimator.output[[j]] = dfj
    }
    df = do.call(rbind, estimator.output)
    df$n = n
    df$p = p
    df$n.test = n.test
    df$dgp = dgp
    df$rho = rho
    df$sim = sim

    out = c(out, list(df))
  }
}
out.df = do.call(rbind, out)

write.csv(out.df, gzfile("simulation.csv.gz"), row.names = FALSE)
