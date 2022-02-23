rm(list = ls())
library(lmtest)
library(grf)
set.seed(123)

n = 2000
p = 15
dgp = "type3"
nreps = 500

# ground truth
data.test = generate_causal_survival_data(50000, p, dgp=dgp, n.mc = 50000)
df = data.frame(cate=data.test$cate, x=data.test$X)
lm1 = coeftest(lm(cate ~ x.1 + x.2, df))
true = lm1[2, 1]
lm2 = coeftest(lm(cate ~ 1, df))
true.ate = lm2[1, 1]

res = replicate(nreps, {
  data = generate_causal_survival_data(n, p, dgp = dgp, n.mc = 1)
  X = data$X
  Y = data$Y
  W = data$W
  D = data$D

  sf = causal_survival_forest(X, Y, W, D, horizon = data$Y.max)
  cate.hat = predict(sf)$pred
  cate.hat.dr = get_scores(sf)

  df = data.frame(cate.hat, cate.hat.dr, x=X)
  lm2 = coeftest(lm(cate.hat ~ x.1 + x.2, df), type="HC3", vcov = sandwich::vcovCL)
  lm3 = coeftest(lm(cate.hat.dr ~ x.1 + x.2, df), type="HC3", vcov = sandwich::vcovCL)
  lm4 = coeftest(lm(cate.hat ~ 1, df), type="HC3", vcov = sandwich::vcovCL)
  lm5 = coeftest(lm(cate.hat.dr ~ 1, df), type="HC3", vcov = sandwich::vcovCL)

  ub = lm2[2, 1] + 1.96*lm2[2, 2]
  lb = lm2[2, 1] - 1.96*lm2[2, 2]
  cov.blp.cate = true > lb & true < ub
  ub = lm3[2, 1] + 1.96*lm3[2, 2]
  lb = lm3[2, 1] - 1.96*lm3[2, 2]
  cov.blp.dr = true > lb & true < ub
  ub = lm4[1, 1] + 1.96*lm4[1, 2]
  lb = lm4[1, 1] - 1.96*lm4[1, 2]
  cov.ate = true.ate > lb & true.ate < ub
  ub = lm5[1, 1] + 1.96*lm5[1, 2]
  lb = lm5[1, 1] - 1.96*lm5[1, 2]
  cov.ate.dr = true.ate > lb & true.ate < ub

  c(blp.cate=lm2[2, 1],
    blp.dr=lm3[2, 1],
    ate=lm4[1, 1],
    ate.dr=lm5[1, 1],
    cov.blp.cate=cov.blp.cate,
    cov.blp.dr=cov.blp.dr,
    cov.ate=cov.ate,
    cov.ate.dr=cov.ate.dr)
})

# Figure 2:
res.cov = round(rowMeans(res), 2)
pdf("blp_simulation.pdf")
breaks = 20
par(mfrow = c(2, 2))
hist(res["blp.cate", ], breaks = breaks, main = "BLP (CATE)", xlab = paste("coverage: ", res.cov["cov.blp.cate"]))
abline(v=true, col = "red", lty = 1)
hist(res["blp.dr", ], breaks = breaks, main = "BLP (DR)", xlab = paste("coverage: ", res.cov["cov.blp.dr"]))
abline(v=true, col = "red", lty = 1)

hist(res["ate", ], breaks = breaks, main = "ATE (CATE)" , xlab = paste("coverage: ", res.cov["cov.ate"]))
abline(v=true.ate, col = "red", lty = 1)
hist(res["ate.dr", ], breaks = breaks, main = "ATE (DR)", xlab = paste("coverage: ", res.cov["cov.ate.dr"]))
abline(v=true.ate, col = "red", lty = 1)
dev.off()
