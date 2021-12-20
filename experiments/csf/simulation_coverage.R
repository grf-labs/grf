rm(list = ls())
library(grf)
set.seed(123)

out = list()
n.sim = 1000
grid = expand.grid(n = c(2000),
                   p = 15,
                   rho = c(0, 0.5),
                   num.trees = c(10000),
                   dgp = c("type1", "type2", "type3", "type4"),
                   stringsAsFactors = FALSE)

for (i in 1:nrow(grid)) {
  print(Sys.time())
  print(paste("grid", i, "of", nrow(grid)))
  print(grid[i, ])
  n = grid$n[i]
  p = grid$p[i]
  dgp = grid$dgp[i]
  rho = grid$rho[i]
  num.trees = grid$num.trees[i]
  X.test = matrix(c(0.2, 0.4, 0.6, 0.8), 4, p)

  data.test = generate_causal_survival_data(n = nrow(X.test), p = p, X = X.test, dgp = dgp, rho = rho, n.mc = 100000)
  cate.true = data.test$cate
  cate.true.prob = data.test$cate.prob
  for (sim in 1:n.sim) {
    print(paste("sim", sim))
    data = generate_causal_survival_data(n = n, p = p, dgp = dgp, rho = rho, n.mc = 1)
    data$Y = round(data$Y, 2)
    forest.W = regression_forest(data$X, data$W, num.trees = 500, ci.group.size = 1)
    W.hat = predict(forest.W)$predictions

    fit = causal_survival_forest(data$X, data$Y, data$W, data$D, W.hat = W.hat, horizon = data$Y.max,
                                 num.trees = num.trees)
    fit.prob = causal_survival_forest(data$X, data$Y, data$W, data$D, W.hat = W.hat,
                                      target = "survival.probability", horizon = data$y0,
                                      num.trees = num.trees)
    pred = predict(fit, X.test, estimate.variance = TRUE)
    pred.prob = predict(fit.prob, X.test, estimate.variance = TRUE)

    lw = pred$predictions - sqrt(pred$variance.estimates) * qnorm(0.975)
    up = pred$predictions + sqrt(pred$variance.estimates) * qnorm(0.975)
    coverage = as.integer(cate.true > lw & cate.true < up)
    width = up - lw
    lw.prob = pred.prob$predictions - sqrt(pred.prob$variance.estimates) * qnorm(0.975)
    up.prob = pred.prob$predictions + sqrt(pred.prob$variance.estimates) * qnorm(0.975)
    coverage.prob = as.integer(cate.true.prob > lw.prob & cate.true.prob < up.prob)
    width.prob = up.prob - lw.prob

    df = data.frame(
      target = c(rep("RMST", 4), rep("survival.probability", 4)),
      coverage = c(coverage, coverage.prob),
      width = c(width, width.prob),
      n = n,
      p = p,
      dgp = dgp,
      rho = rho,
      num.trees = num.trees,
      sim = sim,
      X.test = X.test[, 1]
      )

    out = c(out, list(df))
  }
}
out.df = do.call(rbind, out)

write.csv(out.df, gzfile("coverage.csv.gz"), row.names = FALSE)
