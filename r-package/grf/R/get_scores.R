#' Compute doubly robust scores for a GRF forest object
#'
#' @param forest A grf forest object
#' @param ... Additional arguments
#' @return A vector of scores
#' @export
get_scores <- function(forest, ...) {
  UseMethod("get_scores")
}

#' Compute doubly robust scores for a causal forest.
#'
#' Compute doubly robust (AIPW) scores for average treatment effect estimation
#' or average partial effect estimation with continuous treatment,
#' using a causal forest. Under regularity conditions, the average of the DR.scores
#' is an efficient estimate of the average treatment effect.
#'
#' @param forest A trained causal forest.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) they are obtained via inverse-propensity weighting in the case
#'               of binary treatment or by estimating Var[W | X = x] using a new forest
#'               in the case of a continuous treatment.
#' @param num.trees.for.weights Number of trees used to estimate Var[W | X = x]. Note: this
#'               argument is only used when debiasing.weights = NULL.
#' @param ... Additional arguments (currently ignored).
#'
#' @references Farrell, Max H. "Robust inference on average treatment effects with
#'             possibly more covariates than observations." Journal of Econometrics
#'             189(1), 2015.
#' @references Graham, Bryan S., and Cristine Campos de Xavier Pinto. "Semiparametrically
#'             efficient estimation of the average linear regression function."
#'             Journal of Econometrics 226(1), 2022.
#' @references Hirshberg, David A., and Stefan Wager. "Augmented minimax linear estimation."
#'             The Annals of Statistics 49(6), 2021.
#' @references Robins, James M., and Andrea Rotnitzky. "Semiparametric efficiency in
#'             multivariate regression models with missing data." Journal of the
#'             American Statistical Association 90(429), 1995.
#'
#' @return A vector of scores.
#'
#' @method get_scores causal_forest
#' @export
get_scores.causal_forest <- function(forest,
                                     subset = NULL,
                                     debiasing.weights = NULL,
                                     num.trees.for.weights = 500,
                                     ...) {
  subset <- validate_subset(forest, subset)
  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]

  binary.W <- all(forest$W.orig %in% c(0, 1))

  if (is.null(debiasing.weights)) {
    if (binary.W) {
      debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
    } else {
      # Start by learning debiasing weights if needed.
      # The goal is to estimate the variance of W given X. For binary treatments,
      # we get a good implicit estimator V.hat = e.hat (1 - e.hat), and
      # so this step is not needed. Note that if we use the present CAPE estimator
      # with a binary treatment and set V.hat = e.hat (1 - e.hat), then we recover
      # exactly the AIPW estimator of the CATE.
      clusters <- if (length(forest$clusters) > 0) {
        forest$clusters
      } else {
        1:length(forest$Y.orig)
      }
      variance_forest <- regression_forest(forest$X.orig,
                                           (forest$W.orig - forest$W.hat)^2,
                                           clusters = clusters,
                                           sample.weights = forest$sample.weights,
                                           num.trees = num.trees.for.weights,
                                           ci.group.size = 1,
                                           seed = forest$seed)
      V.hat <- predict(variance_forest)$predictions
      debiasing.weights.all <- (forest$W.orig - forest$W.hat) / V.hat
      debiasing.weights <- debiasing.weights.all[subset]
    }
  } else if (length(debiasing.weights) == length(forest$Y.orig)) {
    debiasing.weights <- debiasing.weights[subset]
  } else if (length(debiasing.weights) != length(subset))  {
    stop("If specified, debiasing.weights must have length n or |subset|.")
  }

  # Form AIPW scores. Note: We are implicitly using the following
  # estimates for the regression surfaces E[Y|X, W=0/1]:
  # Y.hat.0 <- Y.hat - W.hat * tau.hat.pointwise
  # Y.hat.1 <- Y.hat + (1 - W.hat) * tau.hat.pointwise
  Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))

  tau.hat.pointwise + debiasing.weights * Y.residual
}

#' Doubly robust scores for estimating the average conditional local average treatment effect.
#'
#' Given an outcome Y, treatment W and instrument Z, the (conditional) local
#' average treatment effect is tau(x) = Cov[Y, Z | X = x] / Cov[W, Z | X = x].
#' This is the quantity that is estimated with an instrumental forest.
#' It can be intepreted causally in various ways. Given a homogeneity
#' assumption, tau(x) is simply the CATE at x. When W is binary
#' and there are no "defiers", Imbens and Angrist (1994) show that tau(x) can
#' be interpreted as an average treatment effect on compliers. This doubly robust
#' scores provided here are for estimating tau = E[tau(X)].
#'
#' @param forest A trained instrumental forest.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained via the appropriate doubly robust score
#'               construction, e.g., in the case of causal_forests with a binary treatment, they
#'               are obtained via inverse-propensity weighting.
#' @param compliance.score An estimate of the causal effect of Z on W, i.e., Delta(X) = E[W | X, Z = 1]
#'               - E[W | X, Z = 0], which can then be used to produce debiasing.weights. If not provided,
#'               this is estimated via an auxiliary causal forest.
#' @param num.trees.for.weights In some cases (e.g., with causal forests with a continuous
#'               treatment), we need to train auxiliary forests to learn debiasing weights.
#'               This is the number of trees used for this task. Note: this argument is only
#'               used when debiasing.weights = NULL.
#' @param ... Additional arguments (currently ignored).
#'
#' @references Aronow, Peter M., and Allison Carnegie. "Beyond LATE: Estimation of the
#'              average treatment effect with an instrumental variable." Political
#'              Analysis 21(4), 2013.
#' @references Chernozhukov, Victor, Juan Carlos Escanciano, Hidehiko Ichimura,
#'             Whitney K. Newey, and James M. Robins. "Locally robust semiparametric
#'             estimation." Econometrica 90(4), 2022.
#' @references Imbens, Guido W., and Joshua D. Angrist. "Identification and Estimation of
#'             Local Average Treatment Effects." Econometrica 62(2), 1994.
#'
#' @return A vector of scores.
#'
#' @method get_scores instrumental_forest
#' @export
get_scores.instrumental_forest <- function(forest,
                                           subset = NULL,
                                           debiasing.weights = NULL,
                                           compliance.score = NULL,
                                           num.trees.for.weights = 500,
                                           ...) {
  if (!all(forest$Z.orig %in% c(0, 1))) {
    stop(paste(
      "Average conditional local average treatment effect estimation",
      "only implemented for binary instruments."
    ))
  }
  subset <- validate_subset(forest, subset)

  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  Z.orig <- forest$Z.orig[subset]
  Z.hat <- forest$Z.hat[subset]

  if (is.null(debiasing.weights)) {
  # The compliance forest estimates the effect of the "treatment" Z on the "outcome" W.
    if (is.null(compliance.score)) {
      clusters <- if (length(forest$clusters) > 0) {
        forest$clusters
      } else {
        1:length(forest$Y.orig)
      }
      compliance.forest <- causal_forest(X = forest$X.orig,
                                         Y = forest$W.orig,
                                         W = forest$Z.orig,
                                         Y.hat = forest$W.hat,
                                         W.hat = forest$Z.hat,
                                         sample.weights = forest$sample.weights,
                                         clusters = clusters,
                                         num.trees = num.trees.for.weights,
                                         seed = forest$seed)
      compliance.score <- predict(compliance.forest)$predictions
      compliance.score <- compliance.score[subset]
    } else if (length(compliance.score) == length(forest$Y.orig)) {
      compliance.score <- compliance.score[subset]
    } else if (length(compliance.score) != length(subset))  {
      stop("If specified, compliance.score must have length n or |subset|.")
    }
    if (min(Z.hat) <= 0.01 || max(Z.hat) >= 0.99) {
      rng <- range(Z.hat)
      warning(paste0(
        "Estimated instrument propensities take values between ",
        round(rng[1], 3), " and ", round(rng[2], 3),
        " and in particular get very close to 0 or 1. ",
        "Poor overlap may hurt perfmance for average conditional local average ",
        "treatment effect estimation."
      ))
    }
    if (abs(min(compliance.score)) <= 0.01 * sd(W.orig)) {
      warning(paste0(
        "The instrument appears to be weak, with some compliance scores as ",
        "low as ", round(min(compliance.score), 4)
      ))
    }
    debiasing.weights <- (Z.orig - Z.hat) / (Z.hat * (1 - Z.hat)) / compliance.score
  } else if (length(debiasing.weights) == length(forest$Y.orig)) {
    debiasing.weights <- debiasing.weights[subset]
  } else if (length(debiasing.weights) != length(subset))  {
    stop("If specified, debiasing.weights must have length n or |subset|.")
  }

  tau.hat.pointwise <- predict(forest)$predictions[subset]
  Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))

  tau.hat.pointwise + debiasing.weights * Y.residual
}

#' Compute doubly robust scores for a multi arm causal forest.
#'
#' Compute doubly robust (AIPW) scores for average treatment effect estimation
#' using a multi arm causal forest. Under regularity conditions, the average of the DR.scores
#' is an efficient estimate of the average treatment effect.
#'
#' @param forest A trained multi arm causal forest.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param drop If TRUE, coerce the result to the lowest possible dimension. Default is FALSE.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An array of scores for each contrast and outcome.
#'
#' @method get_scores multi_arm_causal_forest
#' @export
get_scores.multi_arm_causal_forest <- function(forest,
                                               subset = NULL,
                                               drop = FALSE,
                                               ...) {
  subset <- validate_subset(forest, subset)
  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset, , drop = FALSE]
  treatment.names <- levels(W.orig)
  observed.treatment <- match(W.orig, treatment.names)
  observed.treatment.idx <- cbind(seq_along(subset), observed.treatment)
  contrast.names <- paste(treatment.names[-1], "-", treatment.names[1])

  if (min(W.hat) <= 0.05 || max(W.hat) >= 0.95) {
    min <- apply(W.hat, 2, min)
    max <- apply(W.hat, 2, max)
    warning(paste0(
      "Estimated treatment propensities take values very close to 0 or 1",
      " meaning some estimates may not be well identified.",
      " In particular, the minimum propensity estimates for each arm is\n",
      paste0(treatment.names, ": ", round(min, 3), collapse = " "),
      "\nand the maximum is\n",
      paste0(treatment.names, ": ", round(max, 3), collapse = " "),
      "."
    ))
  }
  # Fill in a [NxK] IPW matrix with inverse propensity estimates of the observed arm
  # using the subsetting syntax "matrix[index.matrix]".
  IPW <- matrix(0, length(subset), nlevels(W.orig))
  IPW[observed.treatment.idx] <- 1 / W.hat[observed.treatment.idx]
  control <- IPW[, 1] > 0
  IPW[control, -1] <- -1 * IPW[control, 1]

  IPW <- IPW[, -1, drop = FALSE]
  W.hat <- W.hat[, -1, drop = FALSE]
  forest.pp <- predict(forest)

  .get.scores <- function(outcome) {
    Y.orig <- forest$Y.orig[subset, outcome]
    Y.hat <- forest$Y.hat[subset, outcome]
    tau.hat.pointwise <- forest.pp$predictions[subset, , outcome]

    Y.hat.baseline <- Y.hat - rowSums(W.hat * tau.hat.pointwise)
    mu.hat.matrix <- cbind(Y.hat.baseline, Y.hat.baseline + tau.hat.pointwise)
    Y.residual <- Y.orig - mu.hat.matrix[observed.treatment.idx]

    tau.hat.pointwise + IPW * Y.residual
  }

  scores <- lapply(1:NCOL(forest$Y.orig), function(col) .get.scores(col))

  array(unlist(scores),
        dim = c(length(subset), dim(forest.pp$predictions)[-1]),
        dimnames = dimnames(forest.pp$predictions))[, , , drop = drop]
}

#' Compute doubly robust scores for a causal survival forest.
#'
#' For details see section 3.2 in the causal survival forest paper.
#'
#' @param forest A trained causal survival forest.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param num.trees.for.weights Number of trees used to estimate Var[W | X = x]. Note: this
#'  argument is only used in the case of a continuous treatment
#'  (see \code{\link{get_scores.causal_forest}} for details).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A vector of scores.
#'
#' @method get_scores causal_survival_forest
#' @export
get_scores.causal_survival_forest <- function(forest,
                                              subset = NULL,
                                              num.trees.for.weights = 500,
                                              ...) {
  subset <- validate_subset(forest, subset)
  numerator <- forest[["_psi"]]$numerator[subset]
  denominator <- forest[["_psi"]]$denominator[subset]
  cate.hat <- predict(forest)$predictions[subset]
  psi <- numerator - denominator * cate.hat

  binary.W <- all(forest$W.orig %in% c(0, 1))
  if (binary.W) {
    if (min(forest$W.hat[subset]) <= 0.05 || max(forest$W.hat[subset]) >= 0.95) {
      rng <- range(forest$W.hat[subset])
      warning(paste0(
        "Estimated treatment propensities take values very close to 0 or 1.",
        " The estimated propensities are between ",
        round(rng[1], 3), " and ", round(rng[2], 3),
        ", meaning some estimates may not be well identified."
      ))
    }
    W.hat <- forest$W.hat[subset]
    V.hat <- W.hat * (1 - W.hat)
  } else {
    # Estimate Var[W | X = x], as in get_scores.causal_forest.
    clusters <- if (length(forest$clusters) > 0) {
      forest$clusters
    } else {
      1:length(forest$Y.orig)
    }
    variance_forest <- regression_forest(forest$X.orig,
                                         (forest$W.orig - forest$W.hat)^2,
                                         clusters = clusters,
                                         sample.weights = forest$sample.weights,
                                         num.trees = num.trees.for.weights,
                                         ci.group.size = 1,
                                         seed = forest$seed)
    V.hat <- predict(variance_forest)$predictions[subset]
  }

  cate.hat + psi / V.hat
}
