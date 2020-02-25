# dgps.R - convenience script for generating simulation data for grf.
# Usage:
# `source("dgps.R")` exports:
# gen_data() (documented below)


#' Generate comparable simulation data
#'
#' The following function returns one simulated data set from one of the seven following DGPs:
#' "simple", "aw1", "aw2", "aw3", "ai1", "ai2", "kunzel".
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
                     dgp = c("simple", "aw1", "aw2", "aw3", "ai1", "ai2", "kunzel"),
                     ...) {
  .minp <- c(simple=3, aw1=2, aw2=2, aw3=1, ai1=2, ai2=6, kunzel=2)
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
    tau <- 0
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * X[, 1] - 1 + e * tau
    V <- 1
  } else if (dgp == "aw2") {
    # equation (28) of https://arxiv.org/pdf/1510.04342.pdf
    X <- matrix(runif(n * p), n, p)
    zeta1 <- 1 + 1 / (1 + exp(-20 * (X[, 1] - (1 / 3))))
    zeta2 <- 1 + 1 / (1 + exp(-20 * (X[, 2] - (1 / 3))))
    tau <- zeta1 * zeta2
    e <- 0.5
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
    e <- (1 / 4) * (1 + dbeta(X[, 1], 2, 4))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- e * tau
    V <- 1
  } else if (dgp == "ai1") {
    X <- matrix(rnorm(n, p), n, p)
    nu_x <- 0.5 * X[, 1] + X[, 2]
    tau <- 0.25 * X[, 1]
    e <- 0.5
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "ai2") {
    X <- matrix(rnorm(n, p), n, p)
    nu_x <- 0.5 * X[, 1] + 0.5 * X[, 2] + X[, 3] + X[, 4] + X[, 5] + X[, 6]
    tau <- 0.5 * ((X[, 1] > 0) * X[, 1] + (X[, 2] > 0) * X[, 2])
    e <- 0.5
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "kunzel") {
    X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = toeplitz(0.5^seq(0, p - 1)))
    tau <- 8 * (X[, 2] > 0.1)
    beta <- runif(p, -5, 5)
    mu_0 <- X %*% beta + 5 * (X[, 1] > 0.5) + rnorm(n = n)
    mu_1 <- mu_0 + tau + rnorm(n = n)
    e <- 0.01
    W <- rbinom(n = n, size = 1, prob = e)
    m <- W * mu_1 + (1 - W) * mu_0 - (W - e) * tau
    V <- 1
  }
  # Scale and return data (rescale if `m` and `tau` is not constant)
  if (!is.na(sd(m))) {
    m <- m / sd(m) * sigma.m
  }
  if (!is.na(sd(tau))) {
    tau <- tau / sd(tau) * sigma.tau
  }
  V <- V / mean(V) * sigma.noise^2
  Y <- m + (W - e) * tau + sqrt(V) * rnorm(n)
  out <- list(X = X, Y = Y, W = W, tau = tau, m = m, e = e, dgp = dgp)

  out
}
