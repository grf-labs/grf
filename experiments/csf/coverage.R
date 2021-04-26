rm(list = ls())
library(grf)
set.seed(123)

out = list()
n.sim = 200
n.mc = 100000
p = 5
X.test = matrix(c(0.2, 0.4, 0.6, 0.8), 4, p)
grid = expand.grid(n = c(2000),
                   num.trees = c(100000),
                   dgp = c("type1", "type2", "type3", "type4"),
                   stringsAsFactors = FALSE)

for (i in 1:nrow(grid)) {
  print(paste("grid", i))
  n = grid$n[i]
  dgp = grid$dgp[i]
  num.trees = grid$num.trees[i]

  for (sim in 1:n.sim) {
    print(paste("sim", sim))
    data = generate_causal_survival_data(n = n, p = p, dgp = dgp, n.mc = 10)
    data.test = generate_causal_survival_data(n = nrow(X.test), p = p, X = X.test, dgp = dgp, n.mc = n.mc)
    cate.true = data.test$cate
    fit = causal_survival_forest(data$X, data$Y, data$W, data$D,
                                 num.trees = num.trees)
    pred = predict(fit, X.test, estimate.variance = TRUE)

    lw = pred$predictions - sqrt(pred$variance.estimates) * qnorm(0.975)
    up = pred$predictions + sqrt(pred$variance.estimates) * qnorm(0.975)
    coverage = as.integer(cate.true > lw & cate.true < up)
    width = up - lw
    df = data.frame(coverage, width, n, dgp, num.trees, sim, X.test = X.test[, 1])

    out = c(out, list(df))
  }
}
out.df = do.call(rbind, out)

write.csv(out.df, gzfile("coverage.csv.gz"), row.names = FALSE)
