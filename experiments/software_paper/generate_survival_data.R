#' Simulate survival data
#'
#' Survival benchmarks for estimating the conditional expected survival time E[T | X].
#' DGPS "type1" to "type4" are taken from the following source by setting all Wi to 1:
#' \itemize{
#'   \item  "type1": T is drawn from an accelerated failure time model and C from a Cox model (scenario 1 in https://arxiv.org/abs/2001.09887)
#'   \item  "type2": T is drawn from a proportional hazard model and C from a accelerated failure time (scenario 2 in https://arxiv.org/abs/2001.09887)
#'   \item  "type3": T and C are drawn from a Poisson distribution  (scenario 3 in https://arxiv.org/abs/2001.09887)
#'   \item  "type4": T and C are drawn from a Poisson distribution  (scenario 4 in https://arxiv.org/abs/2001.09887)
#' }
#'
#' @param n The number of samples.
#' @param p The number of covariates.
#' @param Y.max The maximum follow-up time (optional).
#' @param X The covariates (optional).
#' @param n.mc The number of monte carlo draws to estimate the treatment effect with. Default is 10000.
#' @param dgp The type of DGP.
#'
#' @return A list with entries:
#'  `X`: the covariates, `Y`: the event times, `D`: the censoring indicator,
#'  `ET`: the expected survival time E[T | X] estimated by monte carlo,
#'  `dgp`: the dgp name, `Y.max`: the maximum follow-up time.
#'
#' @examples
#' \donttest{
#' # Generate data
#' n <- 1000
#' p <- 5
#' data <- generate_survival_data(n, p)
#' # Get true ET on a test set
#' X.test <- matrix(seq(0, 1, length.out = 5), 5, p)
#' data.test <- generate_survival_data(n, p, X = X.test)
#' }
#'
#' @export
generate_survival_data <- function(n, p, Y.max = NULL, X = NULL, n.mc = 10000,
                                   dgp = c("type1", "type2", "type3", "type4")) {
  dgp <- match.arg(dgp)
  if (!is.null(X)) {
    p <- NCOL(X)
    n <- NROW(X)
  }
  # All dgps are the same as those found in grf's `generate_causal_survival_data` with all Wi set to 1.
  if (dgp == "type1") {
    # Type 1 from https://arxiv.org/abs/2001.09887 (Cox PH censor time)
    if (is.null(Y.max)) {
      Y.max <- 1.5
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    W <- rep(1, n)
    I1 <- X[,1 ] < 0.5
    ft <- exp(-1.85 - 0.8 * I1 + 0.7 * sqrt(X[, 2]) + 0.2 * X[, 3] +
                (0.7 - 0.4 * I1 - 0.4 * sqrt(X[, 2])) * W + rnorm(n))
    failure.time <- pmin(ft, Y.max)
    numerator <- -log(runif(n))
    denominator <- exp(-1.75 - 0.5 * sqrt(X[, 2]) + 0.2 * X[, 3]  + (1.15 + 0.5 * I1 - 0.3 * sqrt(X[, 2])) * W)
    censor.time <- (numerator / denominator)^(1/2)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft1 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 2]) + 0.2 * X[i, 3] +
                  0.7 - 0.4 * I1[i] - 0.4 * sqrt(X[i, 2]) + eps)
      ET[i] <- mean(pmin(ft1, Y.max))
    }
  } else if (dgp == "type2") {
    # Type 2 from https://arxiv.org/abs/2001.09887 (Cox PH failure time)
    if (is.null(Y.max)) {
      Y.max <- 2
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    W <- rep(1, n)
    numerator <- -log(runif(n))
    cox.ft <- (numerator / exp(X[,1] + (-0.5 + X[,2]) * W))^2
    failure.time <- pmin(cox.ft, Y.max)
    censor.time <- 3 * runif(n)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft1 <- (numerator / exp(X[i, 1] + (-0.5 + X[i, 2]) * 1))^2
      ET[i] <- mean(pmin(cox.ft1, Y.max))
    }
  } else if (dgp == "type3") {
    # Type 3 from https://arxiv.org/abs/2001.09887 (Poisson)
    if (is.null(Y.max)) {
      Y.max <- 15
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    W <- rep(1, n)
    lambda.failure <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- pmin(rpois(n, lambda = lambda.failure), Y.max)
    lambda.censor <- 12 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3)
    for (i in 1:n) {
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      ET[i] <- mean(pmin(ft1, Y.max))
    }
  } else if (dgp == "type4") {
    # Type 4 from https://arxiv.org/abs/2001.09887 (Poisson)
    if (is.null(Y.max)) {
      Y.max <- 3
    }
    if (is.null(X)) {
      X <- matrix(runif(n * p), n, p)
    }
    W <- rep(1, n)
    lambda.failure <- X[,2] + X[, 3] + pmax(0, X[, 1] - 0.3) * W
    failure.time <- pmin(rpois(n, lambda = lambda.failure), Y.max)
    lambda.censor <- 1 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    lambda.failure.1 <- X[,2] + X[, 3] + pmax(0, X[, 1] - 0.3)
    for (i in 1:n) {
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      ET[i] <- mean(pmin(ft1, Y.max))
    }
  }

  list(X = X, Y = Y, D = D, ET = ET, dgp = dgp, Y.max = Y.max)
}
