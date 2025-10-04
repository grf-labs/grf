set.seed(42)
rm(list = ls())
library(grf)

# Adapted from grf's "generate_causal_survival_data" to benchmark P[T > y0 | X]
generate_data <- function(n, p, Y.max = NULL, y0 = NULL, rho = 0, n.mc = 10000, X = NULL,
                          dgp = c("type1", "type2", "type3", "type4", "type5")) {
  .minp <- c(type1 = 5, type2 = 5, type3 = 5, type4 = 5, type5 = 5)
  dgp <- match.arg(dgp)
  if (p < .minp[dgp]) {
    stop(paste("Selected dgp", dgp, "requires a minimum of", minp, "variables."))
  }
  if (rho !=0 ) {
    if (!("MASS" %in% utils::installed.packages())) {
      stop("`rho != 0` requires the MASS library.")
    }
  }
  if (dgp == "type1") {
    # Type 1 from https://arxiv.org/abs/2001.09887 (Cox PH censor time)
    if (is.null(Y.max)) {
      Y.max <- 1.5
    }
    if (is.null(y0)) {
      y0 <- 0.8 # 90-percentile of Y
    }
    if (is.null(X)) {
      if (rho == 0) {
        X <- matrix(runif(n * p), n, p)
      } else {
        X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
      }
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    I1 <- X[,1 ] < 0.5
    ft <- exp(-1.85 - 0.8 * I1 + 0.7 * sqrt(X[, 2]) + 0.2 * X[, 3] +
                (0.7 - 0.4 * I1 - 0.4 * sqrt(X[, 2])) * W + rnorm(n))
    failure.time <- pmin(ft, Y.max)
    numerator <- -log(runif(n))
    denominator <- exp(-1.75 - 0.5 * sqrt(X[, 2]) + 0.2 * X[, 3] + (1.15 + 0.5 * I1 - 0.3 * sqrt(X[, 2])) * W)
    censor.time <- (numerator / denominator)^(1/2)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    prob <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft0 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 2]) + 0.2 * X[i, 3] + eps)
      ft1 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 2]) + 0.2 * X[i, 3] +
                   0.7 - 0.4 * I1[i] - 0.4 * sqrt(X[i, 2]) + eps)
      prob[i] <- e[i] * mean(ft1 > y0) + (1 - e[i]) * mean(ft0 > y0)
    }
  } else if (dgp == "type2") {
    # Type 2 from https://arxiv.org/abs/2001.09887 (Cox PH failure time)
    if (is.null(Y.max)) {
      Y.max <- 2
    }
    if (is.null(y0)) {
      y0 <- 1.2 # 90-percentile of Y
    }
    if (is.null(X)) {
      if (rho == 0) {
        X <- matrix(runif(n * p), n, p)
      } else {
        X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
      }
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    numerator <- -log(runif(n))
    cox.ft <- (numerator / exp(X[,1] + (-0.5 + X[,2]) * W))^2
    failure.time <- pmin(cox.ft, Y.max)
    censor.time <- 3 * runif(n)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    prob <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft0 <- (numerator / exp(X[i, 1] + (-0.5 + X[i, 2]) * 0))^2
      cox.ft1 <- (numerator / exp(X[i, 1] + (-0.5 + X[i, 2]) * 1))^2
      prob[i] <- e[i] * mean(cox.ft1 > y0) + (1 - e[i]) * mean(cox.ft0 > y0)
    }
  } else if (dgp == "type3") {
    # Type 3 from https://arxiv.org/abs/2001.09887 (Poisson)
    if (is.null(Y.max)) {
      Y.max <- 15
    }
    if (is.null(y0)) {
      y0 <- 10 # 90-percentile of Y
    }
    if (is.null(X)) {
      if (rho == 0) {
        X <- matrix(runif(n * p), n, p)
      } else {
        X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
      }
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    lambda.failure <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- pmin(rpois(n, lambda = lambda.failure), Y.max)
    lambda.censor <- 12 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    prob <- rep(NA, n)
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 6
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      prob[i] <- e[i] * mean(ft1 > y0) + (1 - e[i]) * mean(ft0 > y0)
    }
  } else if (dgp == "type4") {
    # Type 4 from https://arxiv.org/abs/2001.09887 (Poisson)
    if (is.null(Y.max)) {
      Y.max <- 3
    }
    if (is.null(y0)) {
      y0 <- 2 # 90-percentile of Y
    }
    if (is.null(X)) {
      if (rho == 0) {
        X <- matrix(runif(n * p), n, p)
      } else {
        X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
      }
    }
    e <- 1 / ((1 + exp(-X[, 1])) * (1 + exp(-X[, 2])))
    W <- rbinom(n, 1, e)
    lambda.failure <- X[,2] + X[, 3] + pmax(0, X[, 1] - 0.3) * W
    failure.time <- pmin(rpois(n, lambda = lambda.failure), Y.max)
    lambda.censor <- 1 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    prob <- rep(NA, n)
    lambda.failure.0 <- X[,2] + X[, 3]
    lambda.failure.1 <- X[,2] + X[, 3] + pmax(0, X[, 1] - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      prob[i] <- e[i] * mean(ft1 > y0) + (1 - e[i]) * mean(ft0 > y0)
    }
  } else if (dgp == "type5") {
    # Similar to "type2" but censoring generated from an accelerated failure time model.
    if (is.null(Y.max)) {
      Y.max <- 2
    }
    if (is.null(y0)) {
      y0 <- 0.17
    }
    if (is.null(X)) {
      if (rho == 0) {
        X <- matrix(runif(n * p), n, p)
      } else {
        X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
      }
    }
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    numerator <- -log(runif(n))
    cox.ft <- (numerator / exp(X[,1] + (-0.4 + X[,2]) * W))^2
    failure.time <- pmin(cox.ft, Y.max)
    censor.time <- exp(X[1, ] - X[, 3] * W + rnorm(n))
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    prob <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft0 <- (numerator / exp(X[i, 1] + (-0.4 + X[i, 2]) * 0))^2
      cox.ft1 <- (numerator / exp(X[i, 1] + (-0.4 + X[i, 2]) * 1))^2
      prob[i] <- e[i] * mean(cox.ft1 > y0) + (1 - e[i]) * mean(cox.ft0 > y0)
    }
  }

  list(X = X, Y = round(Y, 2), W = W, D = D, dgp = dgp, Y.max = Y.max, y0 = y0,
       true.prob = prob)
}


dgps = c("type1", "type2", "type3", "type4", "type5")
nsim = 250
n = 5000
p = 10

out = list()
for (dgp in dgps) {
  for (sim in 1:nsim) {
    data = generate_data(n = n, p = p, dgp = dgp)
    Y = round(data$Y, 2)
    D = data$D
    X = data$X
    true = data$true.prob

    seed = sim
    sf.exact = survival_forest(X, Y, D, fast.logrank = FALSE, seed = seed, num.trees = 500)
    sf.approx = survival_forest(X, Y, D, fast.logrank = TRUE, seed = seed, num.trees = 500)
    # \hat P[T > h | X]
    pp.exact = predict(sf.exact, prediction.times = "time", failure.times = rep(data$y0, n))$predictions
    pp.approx = predict(sf.approx, prediction.times = "time", failure.times = rep(data$y0, n))$predictions
    # plot(pp.exact, pp.approx); abline(0,1,col='red')
    rmse.exact = sqrt(mean((pp.exact - true)^2))
    rmse.approx = sqrt(mean((pp.approx - true)^2))

    df = data.frame(
      error = rmse.exact - rmse.approx,
      sim = sim,
      dgp = dgp)
    out = c(out, list(df))
  }
}
out.df = do.call(rbind, out)

aggregate(
  list(value = out.df$error),
  by = list(
    name = out.df$dgp
  ),
  FUN = mean
)

# Results
library(ggplot2)

out.df$dgp = factor(out.df$dgp, labels = c("Setting 1", "Setting 2", "Setting 3", "Setting 4", "Setting 5"))
ggplot(out.df, aes(y = error)) +
  geom_boxplot() +
  facet_wrap(~ dgp, nrow = 1) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  ylab(NULL) +
  ggtitle("Difference in RMSE") +
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
  )

ggsave("rmse.pdf", width = 5, height = 3, scale = 1)
write.csv(out.df, "rmse.csv", row.names = FALSE)
