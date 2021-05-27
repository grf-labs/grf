#' Get doubly robust estimates of average treatment effects.
#'
#' In the case of a causal forest with binary treatment, we provide
#' estimates of one of the following:
#' \itemize{
#'   \item The average treatment effect (target.sample = all): E[Y(1) - Y(0)]
#'   \item The average treatment effect on the treated (target.sample = treated): E[Y(1) - Y(0) | Wi = 1]
#'   \item The average treatment effect on the controls (target.sample = control): E[Y(1) - Y(0) | Wi = 0]
#'   \item The overlap-weighted average treatment effect (target.sample = overlap):
#'   E[e(X) (1 - e(X)) (Y(1) - Y(0))] / E[e(X) (1 - e(X)), where e(x) = P[Wi = 1 | Xi = x].
#' }
#' This last estimand is recommended by Li, Morgan, and Zaslavsky (2018)
#' in case of poor overlap (i.e., when the propensities e(x) may be very close
#' to 0 or 1), as it doesn't involve dividing by estimated propensities.
#'
#' In the case of a causal forest with continuous treatment, we provide estimates of the
#' average partial effect, i.e., E[Cov[W, Y | X] / Var[W | X]]. In the case of a binary treatment,
#' the average partial effect matches the average treatment effect. Computing the average partial
#' effect is somewhat more involved, as the relevant doubly robust scores require an estimate
#' of Var[Wi | Xi = x]. By default, we get such estimates by training an auxiliary forest;
#' however, these weights can also be passed manually by specifying debiasing.weights.
#'
#' In the case of instrumental forests with a binary treatment, we provide an estimate
#' of the the Average (Conditional) Local Averate Treatment (ACLATE).
#' Specifically, given an outcome Y, treatment W and instrument Z, the (conditional) local
#' average treatment effect is tau(x) = Cov[Y, Z | X = x] / Cov[W, Z | X = x].
#' This is the quantity that is estimated with an instrumental forest.
#' It can be intepreted causally in various ways. Given a homogeneity
#' assumption, tau(x) is simply the CATE at x. When W is binary
#' and there are no "defiers", Imbens and Angrist (1994) show that tau(x) can
#' be interpreted as an average treatment effect on compliers. This function
#' provides and estimate of tau = E[tau(X)]. See Chernozhukov
#' et al. (2016) for a discussion, and Section 5.2 of Athey and Wager (2021)
#' for an example using forests.
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
#'               Note: Options other than "all" are only currently implemented
#'               for causal forests.
#' @param method Method used for doubly robust inference. Can be either
#'               augmented inverse-propensity weighting (AIPW), or
#'               targeted maximum likelihood estimation (TMLE). Note:
#'               TMLE is currently only implemented for causal forests with
#'               a binary treatment.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained via the appropriate doubly robust score
#'               construction, e.g., in the case of causal_forests with a binary treatment, they
#'               are obtained via inverse-propensity weighting.
#' @param compliance.score Only used with instrumental forests. An estimate of the causal
#'               effect of Z on W, i.e., Delta(X) = E[W | X, Z = 1] - E[W | X, Z = 0],
#'               which can then be used to produce debiasing.weights. If not provided,
#'               this is estimated via an auxiliary causal forest.
#' @param num.trees.for.weights In some cases (e.g., with causal forests with a continuous
#'               treatment), we need to train auxiliary forests to learn debiasing weights.
#'               This is the number of trees used for this task. Note: this argument is only
#'               used when debiasing.weights = NULL.
#'
#' @references Athey, Susan, and Stefan Wager. "Policy Learning With Observational Data."
#'             Econometrica 89.1 (2021): 133-161.
#' @references Chernozhukov, Victor, Juan Carlos Escanciano, Hidehiko Ichimura,
#'             Whitney K. Newey, and James M. Robins. "Locally robust semiparametric
#'             estimation." arXiv preprint arXiv:1608.00033, 2016.
#' @references Imbens, Guido W., and Joshua D. Angrist. "Identification and Estimation of
#'             Local Average Treatment Effects." Econometrica 62(2), 1994.
#' @references Li, Fan, Kari Lock Morgan, and Alan M. Zaslavsky.
#'             "Balancing covariates via propensity score weighting."
#'             Journal of the American Statistical Association 113(521), 2018.
#' @references Robins, James M., and Andrea Rotnitzky. "Semiparametric efficiency in
#'             multivariate regression models with missing data." Journal of the
#'             American Statistical Association 90(429), 1995.
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
#'
#' # Example for causal forests with a continuous treatment.
#' n <- 2000
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(-X[, 2]))) + rnorm(n)
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' tau.forest <- causal_forest(X, Y, W)
#' tau.hat <- predict(tau.forest)
#' average_treatment_effect(tau.forest)
#' average_treatment_effect(tau.forest, subset = X[, 1] > 0)
#' }
#'
#' @return An estimate of the average treatment effect, along with standard error.
#'
#' @export
average_treatment_effect <- function(forest,
                                     target.sample = c("all", "treated", "control", "overlap"),
                                     method = c("AIPW", "TMLE"),
                                     subset = NULL,
                                     debiasing.weights = NULL,
                                     compliance.score = NULL,
                                     num.trees.for.weights = 500) {

  target.sample <- match.arg(target.sample)
  method <- match.arg(method)
  cluster.se <- length(forest$clusters) > 0

  if (method == "TMLE") {
    if (cluster.se) {
      stop("TMLE has not yet been implemented with clustered observations.")
    }
  }

  clusters <- if (cluster.se) {
    forest$clusters
  } else {
    1:NROW(forest$Y.orig)
  }
  observation.weight <- observation_weights(forest)

  subset <- validate_subset(forest, subset)
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  if (!is.null(debiasing.weights)) {
    if (length(debiasing.weights) == NROW(forest$Y.orig)) {
      debiasing.weights <- debiasing.weights[subset]
    } else if (length(debiasing.weights) != length(subset)) {
      stop("If specified, debiasing.weights must be a vector of length n or the subset length.")
    }
  }

  # Add usage guidance for causal forests with a binary treatment, as this
  # is a setting where there are many estimands available

  if ("causal_forest" %in% class(forest) &&
      all(forest$W.orig %in% c(0, 1)) &&
      target.sample != "overlap") {
    if (min(forest$W.hat[subset]) <= 0.05 && max(forest$W.hat[subset]) >= 0.95) {
      rng <- range(forest$W.hat[subset])
      warning(paste0(
        "Estimated treatment propensities take values between ",
        round(rng[1], 3), " and ", round(rng[2], 3),
        " and in particular get very close to 0 and 1. ",
        "In this case, using `target.sample=overlap`, or filtering data as in ",
        "Crump, Hotz, Imbens, and Mitnik (Biometrika, 2009) may be helpful."
      ))
    } else if (min(forest$W.hat[subset]) <= 0.05 && target.sample != "treated") {
      warning(paste0(
        "Estimated treatment propensities go as low as ",
        round(min(forest$W.hat[subset]), 3), " which means that treatment ",
        "effects for some controls may not be well identified. ",
        "In this case, using `target.sample=treated` may be helpful."
      ))
    } else if (max(forest$W.hat[subset]) >= 0.95 && target.sample != "control") {
      warning(paste0(
        "Estimated treatment propensities go as high as ",
        round(max(forest$W.hat[subset]), 3), " which means that treatment ",
        "effects for some treated units may not be well identified. ",
        "In this case, using `target.sample=control` may be helpful."
      ))
    }
  }

  if (method == "AIPW" && target.sample == "all") {
    # This is the most general workflow, that shares codepaths with best linear projection
    # and other average effect estimators.

    .sigma2.hat <- function(DR.scores, tau.hat) {
      correction.clust <- Matrix::sparse.model.matrix(~ factor(subset.clusters) + 0, transpose = TRUE) %*%
        (sweep(as.matrix(DR.scores), 2, tau.hat, "-") * subset.weights)

      Matrix::colSums(correction.clust^2) / sum(subset.weights)^2 *
        nrow(correction.clust) / (nrow(correction.clust) - 1)
    }

    if (any(c("causal_forest", "instrumental_forest", "multi_arm_causal_forest", "causal_survival_forest")
            %in% class(forest))) {
      DR.scores <- get_scores(forest, subset = subset, debiasing.weights = debiasing.weights,
                              compliance.score = compliance.score, num.trees.for.weights = num.trees.for.weights)
    } else {
      stop("Average treatment effects are not implemented for this forest type.")
    }

    if ("multi_arm_causal_forest" %in% class(forest)) {
      tau.hats <- lapply(1:NCOL(forest$Y.orig), function(col) {
        apply(DR.scores[, , col, drop = FALSE], 2, function(dr) weighted.mean(dr, subset.weights))
      })
      sigma2.hats <- lapply(1:NCOL(forest$Y.orig), function(col) {
        .sigma2.hat(DR.scores[, , col], tau.hats[[col]])
      })
      out <- data.frame(
        estimate = unlist(tau.hats),
        std.err = sqrt(unlist(sigma2.hats)),
        contrast = names(unlist(tau.hats)),
        outcome = dimnames(DR.scores)[[3]][rep(1:NCOL(forest$Y.orig), each = dim(DR.scores)[2])],
        stringsAsFactors = FALSE
      )
      return(out) # rownames will be `contrast` when suitable, allowing a convenient `ate["contrast", "estimate"]` access.
    } else {
      tau.hat <- weighted.mean(DR.scores, subset.weights)
      sigma2.hat <- .sigma2.hat(DR.scores, tau.hat)
      return(c(estimate = tau.hat, std.err = sqrt(sigma2.hat)))
    }
  }

  #
  #
  # From here on out, we follow a code path that is specific to causal forests.
  #
  #

  if (!("causal_forest" %in% class(forest))) {
    stop(paste("For any forest type other than causal_forest, the only",
               "implemented option is method=AIPW and target.sample=all"))
  }

  # Only use data selected via subsetting.
  subset.W.orig <- forest$W.orig[subset]
  subset.W.hat <- forest$W.hat[subset]
  subset.Y.orig <- forest$Y.orig[subset]
  subset.Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]

  # Get estimates for the regression surfaces E[Y|X, W=0/1]
  subset.Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
  subset.Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise

  if (target.sample == "overlap") {

    # Address the overlap case separately, as this is a very different estimation problem.
    # The method argument (AIPW vs TMLE) is ignored in this case, as both methods are effectively
    # the same here. Also, note that the overlap-weighted estimator generalizes naturally to the
    # non-binary W case -- see, e.g., Robinson (Econometrica, 1988) -- and so we do not require
    # W to be binary here.

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

  if (!all(forest$W.orig %in% c(0, 1))) {
    stop(paste("With a continuous treatment, only the options target.sample = {all",
               "or overlap} and method=AIPW are implemented."))
  }

  control.idx <- which(subset.W.orig == 0)
  treated.idx <- which(subset.W.orig == 1)

  # Compute naive average effect estimates (notice that this uses OOB)
  if (target.sample == "all") {
    tau.avg.raw <- weighted.mean(tau.hat.pointwise, subset.weights)
    tau.avg.var <-
      sum(subset.weights^2 * (tau.hat.pointwise - tau.avg.raw)^2) /
      sum(subset.weights)^2
  } else if (target.sample == "treated") {
    tau.avg.raw <- weighted.mean(
      tau.hat.pointwise[treated.idx],
      subset.weights[treated.idx]
    )
    tau.avg.var <-
      sum(subset.weights[treated.idx]^2 *
            (tau.hat.pointwise[treated.idx] - tau.avg.raw)^2) /
      sum(subset.weights[treated.idx])^2
  } else if (target.sample == "control") {
    tau.avg.raw <- weighted.mean(
      tau.hat.pointwise[control.idx],
      subset.weights[control.idx]
    )
    tau.avg.var <-
      sum(subset.weights[control.idx]^2 *
            (tau.hat.pointwise[control.idx] - tau.avg.raw)^2) /
      sum(subset.weights[control.idx])^2
  } else {
    stop("Invalid target sample.")
  }

  if (method == "AIPW") {

    # There's no way a user should be able to get here
    if (!(target.sample %in% c("treated", "control"))) {
      stop("Invalid code path")
    }

    # Compute normalized inverse-propensity-type weights gamma
    if (target.sample == "treated") {
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

    dr.correction.all <- subset.W.orig * gamma * (subset.Y.orig - subset.Y.hat.1) -
      (1 - subset.W.orig) * gamma * (subset.Y.orig - subset.Y.hat.0)
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
          B = subset.Y.orig[subset.W.orig == 0] - subset.Y.hat.0[subset.W.orig == 0]
        ))
      eps.tmle.robust.1 <-
        lm(B ~ A + 0, data = data.frame(
          A = 1 / subset.W.hat[subset.W.orig == 1],
          B = subset.Y.orig[subset.W.orig == 1] - subset.Y.hat.1[subset.W.orig == 1]
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
            B = subset.Y.orig[subset.W.orig == 0] - subset.Y.hat.0[subset.W.orig == 0]
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
        ) %*% (subset.Y.orig[subset.W.orig == 1] - subset.Y.hat.1[subset.W.orig == 1])
        s.1 <- sum(delta.1^2) / sum(subset.W.orig == 1) / (sum(subset.W.orig == 1) - 1)
        sigma2.hat <- s.0 + s.1
      } else {
        sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * new.center^2 +
          var(subset.Y.orig[subset.W.orig == 1] - subset.Y.hat.1[subset.W.orig == 1]) / sum(subset.W.orig == 1)
      }
    } else if (target.sample == "control") {
      eps.tmle.robust.1 <-
        lm(B ~ A + 0,
          data = data.frame(
            A = (1 - subset.W.hat[subset.W.orig == 1]) / subset.W.hat[subset.W.orig == 1],
            B = subset.Y.orig[subset.W.orig == 1] - subset.Y.hat.1[subset.W.orig == 1]
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
        ) %*% (subset.Y.orig[subset.W.orig == 0] - subset.Y.hat.0[subset.W.orig == 0])
        s.0 <- sum(delta.0^2) / sum(subset.W.orig == 0) / (sum(subset.W.orig == 0) - 1)
        s.1 <- sandwich::vcovCL(eps.tmle.robust.1, cluster = subset.clusters[subset.W.orig == 1]) *
          new.center^2
        sigma2.hat <- s.0 + s.1
      } else {
        sigma2.hat <- var(subset.Y.orig[subset.W.orig == 0] - subset.Y.hat.0[subset.W.orig == 0]) / sum(subset.W.orig == 0) +
          sandwich::vcovHC(eps.tmle.robust.1) * new.center^2
      }
    } else {
      stop("Invalid target sample.")
    }
  } else {
    stop("Invalid method.")
  }

  tau.avg <- tau.avg.raw + dr.correction
  tau.se <- sqrt(tau.avg.var + sigma2.hat)
  return(c(estimate = unname(tau.avg), std.err = tau.se))
}
