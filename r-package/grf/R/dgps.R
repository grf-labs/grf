# dgps.R - convenience script for generating simulation data for grf.

#' Generate causal forest data
#'
#' The following DGPs are available for benchmarking purposes:
#' \itemize{
#'  \item "simple": tau = max(X1, 0), e = 0.4 + 0.2 * 1{X1 > 0}.
#'  \item "aw1": equation (27) of https://arxiv.org/pdf/1510.04342.pdf
#'  \item "aw2": equation (28) of https://arxiv.org/pdf/1510.04342.pdf
#'  \item "aw3": confounding is from "aw1" and tau is from "aw2"
#'  \item "aw3reverse": Same as aw3, but HTEs anticorrelated with baseline
#'  \item "ai1": "Setup 1" from section 6 of https://arxiv.org/pdf/1504.01132.pdf
#'  \item "ai2": "Setup 2" from section 6 of https://arxiv.org/pdf/1504.01132.pdf
#'  \item "kunzel": "Simulation 1" from A.1 in https://arxiv.org/pdf/1706.03461.pdf
#'  \item "nw1": "Setup A" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'  \item "nw2": "Setup B" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'  \item "nw3": "Setup C" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'  \item "nw4": "Setup D" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
#'}
#'
#' Each DGP is parameterized by
#' X: observables,
#' m: conditional mean of Y,
#' tau: treatment effect,
#' e: propensity scores,
#' V: conditional variance of Y.
#'
#' The following rescaled data is returned
#' m = m / sd(m) * sigma.m,
#' tau = tau / sd(tau) * sigma.tau,
#' V = V / mean(V) * sigma.noise^2,
#' W = rbinom(e),
#' Y = m + (W - e) * tau + sqrt(V) + rnorm(n).
#'
#' @param n The number of observations.
#' @param p The number of covariates (note: the minimum varies by DGP).
#' @param sigma.m The standard deviation of the unconditional mean of Y. Default is 1.
#' @param sigma.tau The standard deviation of the treatment effect. Default  is 0.1.
#' @param sigma.noise The conditional variance of Y. Default is 1.
#' @param dgp The kind of dgp. Default is "simple".
#'
#' @return A list consisting of:
#'  X, Y, W, tau, m, e, dgp.
#'
#' @examples
#' \donttest{
#' # Generate simple benchmark data
#' data <- generate_causal_data(100, 5, dgp = "simple")
#' # Generate data from Wager and Athey (2018)
#' data <- generate_causal_data(100, 5, dgp = "aw1")
#' data2 <- generate_causal_data(100, 5, dgp = "aw2")
#' }
#' @export
generate_causal_data <- function(n, p, sigma.m = 1, sigma.tau = 0.1, sigma.noise = 1,
                                 dgp = c("simple", "aw1", "aw2", "aw3", "aw3reverse",
                                         "ai1", "ai2", "kunzel", "nw1", "nw2", "nw3", "nw4")) {
  # To add an additonal DGP, fill in the template below and add an entry to `dgp` and `.minp`.
  .minp <- c(simple=3, aw1=2, aw2=2, aw3=1, aw3reverse=1,
             ai1=2, ai2=6, kunzel=2, nw1=5, nw2=5, nw3=3, nw4=5)
  dgp <- match.arg(dgp)
  minp <- .minp[dgp]
  if (p < minp) {
    msg <- paste0("Selected dgp ", dgp, " requires a minimum of ", minp, " variables.")
    stop(msg)
  }

  if (dgp == "kunzel") {
    if (!("MASS" %in% utils::installed.packages())) {
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
    X <- matrix(rnorm(n * p), n, p)
    nu_x <- 0.5 * X[, 1] + X[, 2]
    tau <- 0.25 * X[, 1]
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "ai2") {
    X <- matrix(rnorm(n * p), n, p)
    nu_x <- 0.5 * X[, 1] + 0.5 * X[, 2] + X[, 3] + X[, 4] + X[, 5] + X[, 6]
    tau <- 0.5 * ((X[, 1] > 0) * X[, 1] + (X[, 2] > 0) * X[, 2])
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- nu_x + e * tau
    V <- 0.1^2
  } else if (dgp == "kunzel") {
    # "Simulation 1" from A.1 in https://arxiv.org/pdf/1706.03461.pdf
    # Extremely unbalanced treatment assignment, easy treatment effect.
    X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = stats::toeplitz(0.5^seq(0, p - 1)))
    tau <- 8 * (X[, 2] > 0.1)
    beta <- runif(p, -5, 5)
    mu_0 <- X %*% beta + 5 * (X[, 1] > 0.5) + rnorm(n = n)
    mu_1 <- mu_0 + tau + rnorm(n = n)
    e <- rep(0.01, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- c(W * mu_1 + (1 - W) * mu_0 - (W - e) * tau)
    V <- 1
  } else if (dgp == "nw1") {
    # "Setup A" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Difficult nuisance components, easy treatment effect function.
    X <- matrix(runif(n * p), n, p)
    tau <- (X[, 1] + X[, 2]) / 2
    eta <- 0.1
    e <- pmax(eta, pmin(sin(pi * X[, 1] * X[, 2]), 1 - eta))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- sin(pi * X[, 1] * X[, 2]) + 2 * (X[, 3] - 0.5)^2 + X[, 4] + 0.5 * X[, 5] + e * tau
    V <- 1
  } else if (dgp == "nw2") {
    # "Setup B" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Randomized trial
    X <- matrix(rnorm(n * p), n, p)
    tau <- X[,1] + log(1 + exp(X[, 2]))
    e <- rep(0.5, n)
    W <- rbinom(n = n, size = 1, prob = e)
    m <- pmax(0, X[, 1] + X[, 2], X[, 3]) + pmax(0, X[, 4] + X[, 5]) + e * tau
    V <- 1
  } else if (dgp == "nw3") {
    # "Setup C" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Easy propensity score, strong confounding, difficult baseline,
    # constant treatment effect
    X <- matrix(rnorm(n * p), n, p)
    tau <- rep(1, n)
    e <- 1 / (1 + exp(X[, 2] + X[, 3]))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- 2 * log(1 + exp(X[, 1] + X[, 2] + X[, 3])) + e * tau
    V <- 1
  } else if (dgp == "nw4") {
    # "Setup D" from Section 4 of https://arxiv.org/pdf/1712.04912.pdf
    # Unrelated treatment and control arms
    # (No upside to learning them jointly)
    X <- matrix(rnorm(n * p), n, p)
    tau <- pmax(X[, 1] + X[, 2] + X[, 3], 0) - pmax(X[, 4] + X[, 5], 0)
    e <- 1 / (1 + exp(-X[, 1]) + exp(-X[, 2]))
    W <- rbinom(n = n, size = 1, prob = e)
    m <- (pmax(X[, 1] + X[, 2] + X[, 3], 0) + pmax(X[, 4] + X[, 5], 0)) / 2 + e * tau
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

#' Simulate causal survival data
#'
#' The following DGPs are available for benchmarking purposes, T is the failure time
#' and C the censoring time:
#' \itemize{
#'   \item "simple1": T = X1*eps + W, C ~ U(0, 2) where eps ~ Exp(1) and Y.max = 1.
#'   \item  "type1": T is drawn from an accelerated failure time model and C from a Cox model (scenario 1 in https://arxiv.org/abs/2001.09887)
#'   \item  "type2": T is drawn from a proportional hazard model and C from a accelerated failure time (scenario 2 in https://arxiv.org/abs/2001.09887)
#'   \item  "type3": T and C are drawn from a Poisson distribution  (scenario 3 in https://arxiv.org/abs/2001.09887)
#'   \item  "type4": T and C are drawn from a Poisson distribution  (scenario 4 in https://arxiv.org/abs/2001.09887)
#'   \item  "type5": is similar to "type2" but with censoring generated from an accelerated failure time model.
#' }
#' @param n The number of samples.
#' @param p The number of covariates.
#' @param Y.max The maximum follow-up time (optional).
#' @param y0 Query time to estimate P(T(1) > y0 | X) - P(T(0) > y0 | X) (optional).
#' @param X The covariates (optional).
#' @param rho The correlation coefficient of the X's covariance matrix V_{ij} = rho^|i-j|. Default is 0.
#' @param n.mc The number of monte carlo draws to estimate the treatment effect with. Default is 10000.
#' @param dgp The type of DGP.
#'
#' @return A list with entries:
#'  `X`: the covariates, `Y`: the event times, `W`: the treatment indicator, `D`: the censoring indicator,
#'  `cate`: the treatment effect (RMST) estimated by monte carlo, `cate.prob` the difference in survival probability,
#'  `cate.sign`: the true sign of the cate for ITR comparison, `dgp`: the dgp name, `Y.max`: the maximum follow-up time,
#'  `y0`: the query time for difference in survival probability.
#'
#' @examples
#' \donttest{
#' # Generate data
#' n <- 1000
#' p <- 5
#' data <- generate_causal_survival_data(n, p)
#' # Get true CATE on a test set
#' X.test <- matrix(seq(0, 1, length.out = 5), 5, p)
#' cate.test <- generate_causal_survival_data(n, p, X = X.test)$cate
#' }
#'
#' @export
generate_causal_survival_data <- function(n, p, Y.max = NULL, y0 = NULL, X = NULL, rho = 0, n.mc = 10000,
                                          dgp = c("simple1", "type1", "type2", "type3", "type4", "type5")) {
  .minp <- c(simple1 = 1, type1 = 5, type2 = 5, type3 = 5, type4 = 5, type5 = 5)
  dgp <- match.arg(dgp)
  minp <- .minp[dgp]
  if (!is.null(X)) {
    p <- NCOL(X)
    n <- NROW(X)
  }
  if (p < minp) {
    stop(paste("Selected dgp", dgp, "requires a minimum of", minp, "variables."))
  }
  if (rho !=0 ) {
    if (!("MASS" %in% utils::installed.packages())) {
      stop("`rho != 0` requires the MASS library.")
    }
  }

  if (dgp == "simple1") {
    if (is.null(Y.max)) {
      Y.max <- 1
    }
    if (is.null(y0)) {
      y0 <- 0.6
    }
    if (is.null(X)) {
      if (rho == 0) {
        X <- matrix(runif(n * p), n, p)
      } else {
        X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(rho^seq(0, p - 1))))
      }
    }
    W <- rbinom(n, 1, 0.5)
    failure.time <- pmin(rexp(n) * X[, 1] + W, Y.max)
    censor.time <- 2 * runif(n)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    temp <- rexp(n.mc)
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    for (i in 1:n) {
      cate[i] <- mean(pmin(temp * X[i, 1] + 1, Y.max) - pmin(temp * X[i, 1], Y.max))
      cate.prob[i] <- mean(temp * X[i, 1] + 1 > y0) - mean(temp * X[i, 1] > y0)
    }
    cate.sign = rep(1, n)
  } else if (dgp == "type1") {
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
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft0 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 2]) + 0.2 * X[i, 3] + eps)
      ft1 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 2]) + 0.2 * X[i, 3] +
                  0.7 - 0.4 * I1[i] - 0.4 * sqrt(X[i, 2]) + eps)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
      cate.prob[i] <- mean(ft1 > y0) - mean(ft0 > y0)
    }
    cate.sign <- sign(0.7 - 0.4 * I1 - 0.4 * sqrt(X[, 2]))
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
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft0 <- (numerator / exp(X[i, 1] + (-0.5 + X[i, 2]) * 0))^2
      cox.ft1 <- (numerator / exp(X[i, 1] + (-0.5 + X[i, 2]) * 1))^2
      cate[i] <- mean(pmin(cox.ft1, Y.max) - pmin(cox.ft0, Y.max))
      cate.prob[i] <- mean(cox.ft1 > y0) - mean(cox.ft0 > y0)
    }
    cate.sign <- -sign(-0.5 + X[,2]) # Note: negative b/c of Cox model, larger is worse.
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
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 6
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
      cate.prob[i] <- mean(ft1 > y0) - mean(ft0 > y0)
    }
    cate.sign <- sign(sqrt(X[, 1]) - 0.3)
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
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    lambda.failure.0 <- X[,2] + X[, 3]
    lambda.failure.1 <- X[,2] + X[, 3] + pmax(0, X[, 1] - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
      cate.prob[i] <- mean(ft1 > y0) - mean(ft0 > y0)
    }
    cate.sign <- sign(pmax(0, X[, 1] - 0.3))
    # For X1 < 0.3 the cate is zero so both (0, 1) are optimal, and we can ignore this subset.
    cate.sign[X[, 1] < 0.3] <- NA
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
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft0 <- (numerator / exp(X[i, 1] + (-0.4 + X[i, 2]) * 0))^2
      cox.ft1 <- (numerator / exp(X[i, 1] + (-0.4 + X[i, 2]) * 1))^2
      cate[i] <- mean(pmin(cox.ft1, Y.max) - pmin(cox.ft0, Y.max))
      cate.prob[i] <- mean(cox.ft1 > y0) - mean(cox.ft0 > y0)
    }
    cate.sign <- -sign(-0.4 + X[,2]) # Note: negative b/c of Cox model, larger is worse.
  }

  list(X = X, Y = Y, W = W, D = D, cate = cate, cate.prob = cate.prob,
       cate.sign = cate.sign, dgp = dgp, Y.max = Y.max, y0 = y0)
}
