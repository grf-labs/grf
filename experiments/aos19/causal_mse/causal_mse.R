# Table 1 ("C.GRF" column)
# Simulated MSE with the DGPs "aw2" (first rows),  "aw1" (second rows)
# and "aw3" (last rows).
rm(list = ls())
library(grf)

mse.reps <- 60
n <- c(800, 1600)
p <- c(10, 20)
dgp <- c("aw2", "aw1", "aw3")
grid <- expand.grid(n = n, p = p, dgp = dgp, stringsAsFactors = FALSE)

out <- lapply(1:nrow(grid), function(i) {
  n <- grid$n[i]
  p <- grid$p[i]
  dgp <- grid$dgp[i]
  mse <- replicate(mse.reps, {
    data <- generate_causal_data(n = n, p = p, dgp = dgp, sigma.tau = 1)
    data.test <- generate_causal_data(n = 1000, p = p, dgp = dgp, sigma.tau = 1)
    cf <- causal_forest(data$X, data$Y, data$W)
    tau.hat <- predict(cf, data.test$X)$predictions

    mean((tau.hat - data.test$tau)^2)
  })

  data.frame(dgp = dgp, p = p, n = n, C.GRF = mean(mse) * 10)
})

do.call(rbind, out)
