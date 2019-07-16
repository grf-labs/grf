rm(list=ls())
library(grf)

num_sims = 1000
filename = paste0("results",
       Sys.getenv('SLURM_JOB_ID'), "_",
       Sys.getenv('SLURM_LOCALID'), "_",
       Sys.getenv('SLURM_JOB_NAME'), "_",
       as.integer(Sys.time()), ".csv")

for (s in seq(num_sims)) {

  print(paste0("Simulation ", s))

  # Generate data.
  dgp = sample(c("simple"), 1)
  n = sample(c(250, 1000, 5000), 1)
  p = sample(c(10, 20), 1)
  tm = sample(c("earth", "dicekriging"), 1)

  # Create data
  if (dgp == "simple") {
    X = matrix(rnorm(n*p), n, p)
    X.test = matrix(0, 101, p)
    W = rbinom(n, 1, 0.4 + 0.2 * (X[,1] > 0))
    TAU = pmax(X[,1], 0)
    Y = X[,2] + pmin(X[,3], 0) + TAU * W + rnorm(n)
  }

  # Estimate the forest
  cf = causal_forest(X, Y, W, tune.parameters=T, tuning.method=tm)

  # Estimate treatment effect on oob samples
  tau.hat.oob = predict(cf)$predictions

  # Compute mse
  mse.oob = mean((TAU - tau.hat.oob)^2)

  # Save results
  res = rbind(c(cf$tuning.output$params, tuning.method=tm, mse.oob=mse.oob))
  write.table(res, file=filename, col.names=s == 0, row.names=F, append=T)

}
