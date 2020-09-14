# dgps.R - convenience script for generating simulation data for grf.
# Usage:
# `source("dgps.R")` exports:
# gen_data() (documented below)


#' Generate comparable simulation data
#'
#' The following function returns one simulated data set from one of the eleven following DGPs:
#' "simple", "aw1", "aw2", "aw3", "ai1", "ai2", "kunzel", "nw1", "nw2", "nw3", "nw4"
#'
#' Each DGP is parameterized by
#' x: observables
#' m: conditional mean of Y
#' tau: treatment effect
#' e: propensity scores
#' V: conditional variance of Y
#'
#' The following rescaled data is returned
#' m = m / sd(m) * sigma.m
#' tau = tau / sd(tau) * sigma.tau
#' V = V / mean(V) * sigma.noise^2
#' W = rbinom(e)
#' Y = m + (W - e) * tau + sqrt(V) + rnorm(n)
#'
#' @param n The number of observations
#' @param p The number of covariates (note: the minimum varies by DGP)
#' @param sigma.m The standard deviation of the unconditional mean of Y. Default=1.
#' @param sigma.tau The standard deviation of the treatment effect. Default=0.1.
#' @param sigma.noise The conditional variance of Y. Default=1.
#' @param dgp The kind of dgp. Default="simple".
#'
#' @return A list consisting of:
#'  X, Y, W, tau, m, e, dgp.
#'
#' @examples
#' gen_data(100, 5, sigma.m = 1, sigma.tau = 0.1, sigma.noise = 1, dgp = "simple")
#' @note:
#' To add an additonal DGP, fill in the template below and add an entry to `dgp` and `.minp`.
#'
gen_data <- function(n, p, sigma.m = 1, sigma.tau = 0.1, sigma.noise = 1,
                     dgp = c("simple", "aw1", "aw2", "aw3", "ai1", "ai2", "kunzel",
                             "nw1", "nw2", "nw3", "nw4"),
                     ...) {
  .minp <- c(simple=3, aw1=2, aw2=2, aw3=1, ai1=2, ai2=6, kunzel=2,
             nw1=5, nw2=5, nw3=3, nw4=5)
  dgp <- match.arg(dgp)
  minp <- .minp[dgp]
  if (p < minp) {
    msg <- paste0("Selected dgp ", dgp, " requires a minimum of ", minp, " variables.")
    stop(msg)
  }

  if (dgp == "kunzel") {
    if (!("MASS" %in% installed.packages())) {
      msg <- paste0("Selected dgp ", dgp, " requires the MASS library.")
      stop(msg)
    }
  }

  # Create data
  if (dgp == "simple") {
    X <- matrix(rnorm(n * p), n, p)
    tau <- pmax(X[, 1], 0)
    e <- 0.4 + 0.2 * (X[, 1] > 0)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- X[, 2] + pmin(X[, 3], 0) + e * tau
    V <- 1
  } else if (dgp == "aw1") {
    # equation (27) of https://arxiv.org/pdf/1510.04342.pdf
    X <- matrix(runif(n * p, min = 0, max = 1), n, p)
    tau <- rep(0, n)  # Treatment effect is zero
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))  # Confounding
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * X[, 1] - 1 + e * tau
    V <- 1
  } else if (dgp == "aw2") {
    # equation (28) of https://arxiv.org/pdf/1510.04342.pdf
    X <- matrix(runif(n * p), n, p)
    zeta1 <- 1 + 1 / (1 + exp(-20 * (X[, 1] - (1 / 3))))
    zeta2 <- 1 + 1 / (1 + exp(-20 * (X[, 2] - (1 / 3))))
    tau <- zeta1 * zeta2
    e <- rep(0.5, n)  # Randomized trial (no confounding)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- e * tau
    V <- 1
  } else if (dgp == "aw3") {
    # section 6.2 in https://arxiv.org/pdf/1610.01271.pdf
    # (confounding from aw1, tau from aw2)
    X <- matrix(runif(n * p), n, p)
    zeta1 <- 1 + 1 / (1 + exp(-20 * (X[, 1] - (1 / 3))))
    zeta2 <- 1 + 1 / (1 + exp(-20 * (X[, 2] - (1 / 3))))
    tau <- zeta1 * zeta2
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))  # Confounding
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * X[, 1] - 1 + e * tau
    V <- 1
  } else if (dgp == "aw3reverse") {
    # Same as aw3, but HTEs anticorrelated with baseline
    X <- matrix(runif(n * p), n, p)
    zeta1 <- 1 + 1 / (1 + exp(20 * (X[, 1] - (1 / 3))))
    zeta2 <- 1 + 1 / (1 + exp(20 * (X[, 2] - (1 / 3))))
    tau <- zeta1 * zeta2
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))  # Confounding
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * X[, 1] - 1 + e * tau
    V <- 1
  } else if (dgp == "ai1") {
    X <- matrix(rnorm(n, p), n, p)
    nu_x <- 0.5 * X[, 1] + X[, 2]
    tau <- 0.25 * X[, 1]
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "ai2") {
    X <- matrix(rnorm(n, p), n, p)
    nu_x <- 0.5 * X[, 1] + 0.5 * X[, 2] + X[, 3] + X[, 4] + X[, 5] + X[, 6]
    tau <- 0.5 * ((X[, 1] > 0) * X[, 1] + (X[, 2] > 0) * X[, 2])
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "kunzel") {
    # "Simulation 1" from A.1 in https://arxiv.org/pdf/1706.03461.pdf
    # Extremely unbalanced treatment assignment, easy treatment effect.
    X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = toeplitz(0.5^seq(0, p - 1)))
    tau <- 8 * (X[, 2] > 0.1)
    beta <- runif(p, -5, 5)
    mu_0 <- X %*% beta + 5 * (X[, 1] > 0.5) + rnorm(n = n)
    mu_1 <- mu_0 + tau + rnorm(n = n)
    e <- rep(0.01, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- W * mu_1 + (1 - W) * mu_0 - (W - e) * tau
    V <- 1
  } else if (dgp == "nw1") {
    # "Setup A" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Difficult nuisance components, easy treatment effect function.
    X <- matrix(runif(n * p), n, p)
    tau <- (X[, 1] + X[, 2]) / 2
    eta <- 0.1
    e <- pmax(eta, pmin(sin(pi * X[, 1] * X[, 2]), 1 - eta))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- sin(pi * X[, 1] * X[, 2]) + 2 * (X[, 3] - 0.5)^2 + X[, 4] + 0.5 * X[, 5]
    V <- 1
  } else if (dgp == "nw2") {
    # "Setup B" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Randomized trial
    X <- matrix(rnorm(n * p), n, p)
    tau <- X[,1] + log(1 + exp(X[, 2]))
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- pmax(0, X[, 1] + X[, 2], X[, 3]) + pmax(0, X[, 4] + X[, 5])
    V <- 1
  } else if (dgp == "nw3") {
    # "Setup C" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Easy propensity score, strong confounding, difficult baseline,
    # constant treatment effect
    X <- matrix(rnorm(n * p), n, p)
    tau <- rep(1, n)
    e <- 1 / (1 + exp(X[, 2] + X[, 3]))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * log(1 + exp(X[, 1] + X[, 2] + X[, 3]))
    V <- 1
  } else if (dgp == "nw4") {
    # "Setup D" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Unrelated treatment and control arms
    # (No upside to learning them jointly)
    X <- matrix(rnorm(n * p), n, p)
    tau <- pmax(X[, 1] + X[, 2] + X[, 3], 0) - pmax(X[, 4] + X[, 5], 0)
    e <- 1 / (1 + exp(-X[, 1]) + exp(-X[, 2]))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- (pmax(X[, 1] + X[, 2] + X[, 3], 0) + pmax(X[, 4] + X[, 5], 0)) / 2
    V <- 1
  }

  # Scale and return data (rescale if `m` and `tau` is not constant, the NA check is for when n=1)
  if (!is.na(sd(m)) & !(sd(m) == 0)) {
    m <- m / sd(m) * sigma.m
  }
  if (!is.na(sd(tau)) & !(sd(tau) == 0)) {
    tau <- tau / sd(tau) * sigma.tau
  }
  V <- V / mean(V) * sigma.noise^2
  Y <- m + (W - e) * tau + sqrt(V) * rnorm(n)
  out <- list(X = X, Y = Y, W = W, tau = tau, m = m, e = e, dgp = dgp)

  out
}
