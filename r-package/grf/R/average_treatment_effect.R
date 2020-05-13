#' Estimate average treatment effects using a causal forest
#'
#' Gets estimates of one of the following.
#' \itemize{
#'   \item The (conditional) average treatment effect (target.sample = all):
#'   sum_{i = 1}^n E[Y(1) - Y(0) | X = Xi] / n
#'   \item The (conditional) average treatment effect on the treated (target.sample = treated):
#'   sum_{Wi = 1} E[Y(1) - Y(0) | X = Xi] / |{i : Wi = 1}|
#'   \item The (conditional) average treatment effect on the controls (target.sample = control):
#'   sum_{Wi = 0} E[Y(1) - Y(0) | X = Xi] / |{i : Wi = 0}|
#'   \item The overlap-weighted (conditional) average treatment effect
#'   sum_{i = 1}^n e(Xi) (1 - e(Xi)) E[Y(1) - Y(0) | X = Xi] / sum_{i = 1}^n e(Xi) (1 - e(Xi)),
#'   where e(x) = P[Wi = 1 | Xi = x].
#' }
#' This last estimand is recommended by Li, Morgan, and Zaslavsky (JASA, 2017)
#' in case of poor overlap (i.e., when the propensities e(x) may be very close
#' to 0 or 1), as it doesn't involve dividing by estimated propensities.
#'
#' If clusters are specified, then each unit gets equal weight by default. For
#' example, if there are 10 clusters with 1 unit each and per-cluster ATE = 1,
#' and there are 10 clusters with 19 units each and per-cluster ATE = 0, then
#' the overall ATE is 0.05 (additional sample.weights allow for custom
#' weighting). If equalize.cluster.weights = TRUE each cluster gets equal weight
#' and the overall ATE is 0.5.
#'
#' @param forest The trained forest.
#' @param target.sample Which sample to aggregate treatment effects over.
#' @param method Method used for doubly robust inference. Can be either
#'               augmented inverse-propensity weighting (AIPW), or
#'               targeted maximum likelihood estimation (TMLE).
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#'
#' @examples
#' \donttest{
#' # Train a causal forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' c.forest <- causal_forest(X, Y, W)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' c.pred <- predict(c.forest, X.test)
#' # Estimate the conditional average treatment effect on the full sample (CATE).
#' average_treatment_effect(c.forest, target.sample = "all")
#'
#' # Estimate the conditional average treatment effect on the treated sample (CATT).
#' # We don't expect much difference between the CATE and the CATT in this example,
#' # since treatment assignment was randomized.
#' average_treatment_effect(c.forest, target.sample = "treated")
#'
#' # Estimate the conditional average treatment effect on samples with positive X[,1].
#' average_treatment_effect(c.forest, target.sample = "all", subset = X[, 1] > 0)
#' }
#'
#' @return An estimate of the average treatment effect, along with standard error.
#'
#' @importFrom stats coef lm predict var weighted.mean
#' @export
average_treatment_effect <- function(forest,
                                     target.sample = c("all", "treated", "control", "overlap"),
                                     method = c("AIPW", "TMLE"),
                                     subset = NULL) {
  target.sample <- match.arg(target.sample)
  method <- match.arg(method)
  cluster.se <- length(forest$clusters) > 0

  if (!("causal_forest" %in% class(forest))) {
    stop("Average effect estimation only implemented for causal_forest")
  }

  if (cluster.se && method == "TMLE") {
    stop("TMLE has not yet been implemented with clustered observations.")
  }

  if (is.null(subset)) {
    subset <- 1:length(forest$Y.hat)
  }

  if (class(subset) == "logical" && length(subset) == length(forest$Y.hat)) {
    subset <- which(subset)
  }

  if (!all(subset %in% 1:length(forest$Y.hat))) {
    stop(paste(
      "If specified, subset must be a vector contained in 1:n,",
      "or a boolean vector of length n."
    ))
  }

  clusters <- if (cluster.se) {
    forest$clusters
  } else {
    1:length(forest$Y.orig)
  }
  observation.weight <- observation_weights(forest)

  # Only use data selected via subsetting.
  subset.W.orig <- forest$W.orig[subset]
  subset.W.hat <- forest$W.hat[subset]
  subset.Y.orig <- forest$Y.orig[subset]
  subset.Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  # Address the overlap case separately, as this is a very different estimation problem.
  # The method argument (AIPW vs TMLE) is ignored in this case, as both methods are effectively
  # the same here. Also, note that the overlap-weighted estimator generalizes naturally to the
  # non-binary W case -- see, e.g., Robinson (Econometrica, 1988) -- and so we do not require
  # W to be binary here.

  if (target.sample == "overlap") {
    W.residual <- subset.W.orig - subset.W.hat
    Y.residual <- subset.Y.orig - subset.Y.hat
    tau.ols <- lm(Y.residual ~ W.residual, weights = subset.weights)
    tau.est <- coef(summary(tau.ols))[2, 1]

    if (cluster.se) {
      tau.se <- sqrt(sandwich::vcovCL(tau.ols, cluster = subset.clusters)[2, 2])
    } else {
      tau.se <- sqrt(sandwich::vcovHC(tau.ols)[2, 2])
    }

    return(c(estimate = tau.est, std.err = tau.se))
  }

  if (!all(subset.W.orig %in% c(0, 1))) {
    stop(paste(
      "Average treatment effect estimation only implemented for binary treatment.",
      "See `average_partial_effect` for continuous W."
    ))
  }

  if (min(subset.W.hat) <= 0.01 && max(subset.W.hat) >= 0.99) {
    rng <- range(subset.W.hat)
    warning(paste0(
      "Estimated treatment propensities take values between ",
      round(rng[1], 3), " and ", round(rng[2], 3),
      " and in particular get very close to 0 and 1. ",
      "In this case, using `target.sample=overlap`, or filtering data as in ",
      "Crump, Hotz, Imbens, and Mitnik (Biometrika, 2009) may be helpful."
    ))
  } else if (min(subset.W.hat) <= 0.01 && target.sample != "treated") {
    warning(paste0(
      "Estimated treatment propensities go as low as ",
      round(min(subset.W.hat), 3), " which means that treatment ",
      "effects for some controls may not be well identified. ",
      "In this case, using `target.sample=treated` may be helpful."
    ))
  } else if (max(subset.W.hat) >= 0.99 && target.sample != "control") {
    warning(paste0(
      "Estimated treatment propensities go as high as ",
      round(max(subset.W.hat), 3), " which means that treatment ",
      "effects for some treated units may not be well identified. ",
      "In this case, using `target.sample=control` may be helpful."
    ))
  }

  control.idx <- which(subset.W.orig == 0)
  treated.idx <- which(subset.W.orig == 1)

  # Compute naive average effect estimates (notice that this uses OOB)
  if (target.sample == "all") {
    tau.avg.raw <- weighted.mean(tau.hat.pointwise, subset.weights)
  } else if (target.sample == "treated") {
    tau.avg.raw <- weighted.mean(
      tau.hat.pointwise[treated.idx],
      subset.weights[treated.idx]
    )
  } else if (target.sample == "control") {
    tau.avg.raw <- weighted.mean(
      tau.hat.pointwise[control.idx],
      subset.weights[control.idx]
    )
  } else {
    stop("Invalid target sample.")
  }

  # Get estimates for the regress surfaces E[Y|X, W=0/1]
  Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
  Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise

  if (method == "TMLE") {
    loaded <- requireNamespace("sandwich", quietly = TRUE)
    if (!loaded) {
      warning("To use TMLE, please install the package `sandwich`. Using AIPW instead.")
      method <- "AIPW"
    }
  }

  # Now apply a doubly robust correction
  if (method == "AIPW") {

    # Compute normalized inverse-propensity-type weights gamma
    if (target.sample == "all") {
      gamma.control.raw <- 1 / (1 - subset.W.hat[control.idx])
      gamma.treated.raw <- 1 / subset.W.hat[treated.idx]
    } else if (target.sample == "treated") {
      gamma.control.raw <- subset.W.hat[control.idx] / (1 - subset.W.hat[control.idx])
      gamma.treated.raw <- rep(1, length(treated.idx))
    } else if (target.sample == "control") {
      gamma.control.raw <- rep(1, length(control.idx))
      gamma.treated.raw <- (1 - subset.W.hat[treated.idx]) / subset.W.hat[treated.idx]
    } else {
      stop("Invalid target sample.")
    }

    gamma <- rep(0, length(subset.W.orig))
    gamma[control.idx] <- gamma.control.raw /
      sum(subset.weights[control.idx] * gamma.control.raw) *
      sum(subset.weights)
    gamma[treated.idx] <- gamma.treated.raw /
      sum(subset.weights[treated.idx] * gamma.treated.raw) *
      sum(subset.weights)

    dr.correction.all <- subset.W.orig * gamma * (subset.Y.orig - Y.hat.1) -
      (1 - subset.W.orig) * gamma * (subset.Y.orig - Y.hat.0)
    dr.correction <- weighted.mean(dr.correction.all, subset.weights)

    if (cluster.se) {
      correction.clust <- Matrix::sparse.model.matrix(
        ~ factor(subset.clusters) + 0,
        transpose = TRUE
      ) %*% (dr.correction.all * subset.weights)
      sigma2.hat <- sum(correction.clust^2) / sum(subset.weights)^2 *
        length(correction.clust) / (length(correction.clust) - 1)
    } else {
      sigma2.hat <- mean(dr.correction.all^2) / (length(dr.correction.all) - 1)
    }
  } else if (method == "TMLE") {
    if (target.sample == "all") {
      eps.tmle.robust.0 <-
        lm(B ~ A + 0, data = data.frame(
          A = 1 / (1 - subset.W.hat[subset.W.orig == 0]),
          B = subset.Y.orig[subset.W.orig == 0] - Y.hat.0[subset.W.orig == 0]
        ))
      eps.tmle.robust.1 <-
        lm(B ~ A + 0, data = data.frame(
          A = 1 / subset.W.hat[subset.W.orig == 1],
          B = subset.Y.orig[subset.W.orig == 1] - Y.hat.1[subset.W.orig == 1]
        ))
      delta.tmle.robust.0 <- predict(eps.tmle.robust.0, newdata = data.frame(A = mean(1 / (1 - subset.W.hat))))
      delta.tmle.robust.1 <- predict(eps.tmle.robust.1, newdata = data.frame(A = mean(1 / subset.W.hat)))
      dr.correction <- delta.tmle.robust.1 - delta.tmle.robust.0
      # use robust SE
      if (cluster.se) {
        sigma2.hat <- sandwich::vcovCL(eps.tmle.robust.0, cluster = subset.clusters[subset.W.orig == 0]) *
          mean(1 / (1 - subset.W.hat))^2 +
          sandwich::vcovCL(eps.tmle.robust.1, cluster = subset.clusters[subset.W.orig == 1]) *
            mean(1 / subset.W.hat)^2
      } else {
        sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * mean(1 / (1 - subset.W.hat))^2 +
          sandwich::vcovHC(eps.tmle.robust.1) * mean(1 / subset.W.hat)^2
      }
    } else if (target.sample == "treated") {
      eps.tmle.robust.0 <-
        lm(B ~ A + 0,
          data = data.frame(
            A = subset.W.hat[subset.W.orig == 0] / (1 - subset.W.hat[subset.W.orig == 0]),
            B = subset.Y.orig[subset.W.orig == 0] - Y.hat.0[subset.W.orig == 0]
          )
        )
      new.center <- mean(subset.W.hat[subset.W.orig == 1] / (1 - subset.W.hat[subset.W.orig == 1]))
      delta.tmle.robust.0 <- predict(eps.tmle.robust.0,
        newdata = data.frame(A = new.center)
      )
      dr.correction <- -delta.tmle.robust.0
      if (cluster.se) {
        s.0 <- sandwich::vcovCL(eps.tmle.robust.0, cluster = subset.clusters[subset.W.orig == 0]) *
          new.center^2
        delta.1 <- Matrix::sparse.model.matrix(
          ~ factor(subset.clusters[subset.W.orig == 1]) + 0,
          transpose = TRUE
        ) %*% (subset.Y.orig[subset.W.orig == 1] - Y.hat.1[subset.W.orig == 1])
        s.1 <- sum(delta.1^2) / sum(subset.W.orig == 1) / (sum(subset.W.orig == 1) - 1)
        sigma2.hat <- s.0 + s.1
      } else {
        sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * new.center^2 +
          var(subset.Y.orig[subset.W.orig == 1] - Y.hat.1[subset.W.orig == 1]) / sum(subset.W.orig == 1)
      }
    } else if (target.sample == "control") {
      eps.tmle.robust.1 <-
        lm(B ~ A + 0,
          data = data.frame(
            A = (1 - subset.W.hat[subset.W.orig == 1]) / subset.W.hat[subset.W.orig == 1],
            B = subset.Y.orig[subset.W.orig == 1] - Y.hat.1[subset.W.orig == 1]
          )
        )
      new.center <- mean((1 - subset.W.hat[subset.W.orig == 0]) / subset.W.hat[subset.W.orig == 0])
      delta.tmle.robust.1 <- predict(eps.tmle.robust.1,
        newdata = data.frame(A = new.center)
      )
      dr.correction <- delta.tmle.robust.1
      if (cluster.se) {
        delta.0 <- Matrix::sparse.model.matrix(
          ~ factor(subset.clusters[subset.W.orig == 0]) + 0,
          transpose = TRUE
        ) %*% (subset.Y.orig[subset.W.orig == 0] - Y.hat.0[subset.W.orig == 0])
        s.0 <- sum(delta.0^2) / sum(subset.W.orig == 0) / (sum(subset.W.orig == 0) - 1)
        s.1 <- sandwich::vcovCL(eps.tmle.robust.1, cluster = subset.clusters[subset.W.orig == 1]) *
          new.center^2
        sigma2.hat <- s.0 + s.1
      } else {
        sigma2.hat <- var(subset.Y.orig[subset.W.orig == 0] - Y.hat.0[subset.W.orig == 0]) / sum(subset.W.orig == 0) +
          sandwich::vcovHC(eps.tmle.robust.1) * new.center^2
      }
    } else {
      stop("Invalid target sample.")
    }
  } else {
    stop("Invalid method.")
  }

  tau.avg <- tau.avg.raw + dr.correction
  tau.se <- sqrt(sigma2.hat)
  return(c(estimate = tau.avg, std.err = tau.se))
}
