#' Omnibus evaluation of the quality of the random forest estimates via calibration.
#'
#' Test calibration of the forest. Computes the best linear fit of the target
#' estimand using the forest prediction (on held-out data) as well as the mean
#' forest prediction as the sole two regressors. A coefficient of 1 for
#' `mean.forest.prediction` suggests that the mean forest prediction is correct,
#' whereas a coefficient of 1 for `differential.forest.prediction` additionally suggests
#' that the forest has captured heterogeneity in the underlying signal.
#' The p-value of the `differential.forest.prediction` coefficient
#' also acts as an omnibus test for the presence of heterogeneity: If the coefficient
#' is significantly greater than 0, then we can reject the null of
#' no heterogeneity.
#'
#' @param forest The trained forest.
#' @return A heteroskedasticity-consistent test of calibration.
#' @references Chernozhukov, Victor, Mert Demirer, Esther Duflo, and Ivan Fernandez-Val.
#'             "Generic Machine Learning Inference on Heterogenous Treatment Effects in
#'             Randomized Experiments." arXiv preprint arXiv:1712.04802 (2017).
#'
#' @examples
#' \dontrun{
#' n <- 800
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' forest <- causal_forest(X, Y, W)
#' test_calibration(forest)
#' }
#'
#' @export
test_calibration <- function(forest) {
  observation.weight <- observation_weights(forest)
  clusters <- if (length(forest$clusters) > 0) {
    forest$clusters
  } else {
    1:length(observation.weight)
  }
  if ("regression_forest" %in% class(forest)) {
    preds <- predict(forest)$predictions
    mean.pred <- weighted.mean(preds, observation.weight)
    DF <- data.frame(
      target = unname(forest$Y.orig),
      mean.forest.prediction = mean.pred,
      differential.forest.prediction = preds - mean.pred
    )
  } else if ("causal_forest" %in% class(forest)) {
    preds <- predict(forest)$predictions
    mean.pred <- weighted.mean(preds, observation.weight)
    DF <- data.frame(
      target = unname(forest$Y.orig - forest$Y.hat),
      mean.forest.prediction = unname(forest$W.orig - forest$W.hat) * mean.pred,
      differential.forest.prediction = unname(forest$W.orig - forest$W.hat) *
        (preds - mean.pred)
    )
  } else {
    stop("Calibration check not supported for this type of forest.")
  }

  best.linear.predictor <-
    lm(target ~ mean.forest.prediction + differential.forest.prediction + 0,
      weights = observation.weight,
      data = DF
    )
  blp.summary <- lmtest::coeftest(best.linear.predictor,
    vcov = sandwich::vcovCL,
    type = "HC3",
    cluster = clusters
  )
  attr(blp.summary, "method") <-
    paste("Best linear fit using forest predictions (on held-out data)",
      "as well as the mean forest prediction as regressors, along",
      "with one-sided heteroskedasticity-robust (HC3) SEs",
      sep = "\n"
    )
  # convert to one-sided p-values
  dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
  blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[, 4] / 2, blp.summary[, 4] / 2)
  blp.summary
}



#' Estimate the best linear projection of a conditional average treatment effect
#' using a causal forest.
#'
#' Let tau(Xi) = E[Y(1) - Y(0) | X = Xi] be the CATE, and Ai be a vector of user-provided
#' covariates. This function provides a (doubly robust) fit to the linear model
#'
#' tau(Xi) ~ beta_0 + Ai * beta
#'
#' Procedurally, we do so by regressing doubly robust scores derived from the causal
#' forest against the Ai. Note the covariates Ai may consist of a subset of the Xi,
#' or they may be distinct The case of the null model tau(Xi) ~ beta_0 is equivalent
#' to fitting an average treatment effect via AIPW.
#'
#' In the event the treatment is continuous the inverse-propensity weight component of the
#' double robust scores are replaced with a component based on a forest based
#' estimate of Var[Wi | Xi = x]. These weights can also be passed manually by specifying
#' debiasing.weights.
#'
#' @param forest The trained forest.
#' @param A The covariates we want to project the CATE onto.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n of debiasing weights. If NULL (default)
#'                          and the treatment is binary, then inverse-propensity weighting is used,
#'                          otherwise, if the treatment is continuous, these are estimated by a variance
#'                          forest.
#' @param num.trees.for.variance Number of trees used to estimate Var[Wi | Xi = x]. Default is 500.
#'                               (only applies with continuous treatment and debiasing.weights = NULL)
#'
#' @references Chernozhukov, Victor, and Vira Semenova. "Simultaneous inference for
#'             Best Linear Predictor of the Conditional Average Treatment Effect and
#'             other structural functions." arXiv preprint arXiv:1702.06240 (2017).
#'
#' @examples
#' \dontrun{
#' n <- 800
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' forest <- causal_forest(X, Y, W)
#' best_linear_projection(forest, X[,1:2])
#' }
#'
#' @return An estimate of the best linear projection, along with coefficient standard errors.
#'
#' @importFrom stats lm
#'
#' @export
best_linear_projection <- function(forest,
                                   A = NULL,
                                   subset = NULL,
                                   debiasing.weights = NULL,
                                   num.trees.for.variance = 500) {
  if (!("causal_forest" %in% class(forest))) {
    stop("`best_linear_projection` is only implemented for `causal_forest`")
  }

  if (is.null(subset)) {
    subset <- 1:length(forest$Y.hat)
  }

  if (class(subset) == "logical" & length(subset) == length(forest$Y.hat)) {
    subset <- which(subset)
  }

  if (!all(subset %in% 1:length(forest$Y.hat))) {
    stop(paste(
      "If specified, subset must be a vector contained in 1:n,",
      "or a boolean vector of length n."
    ))
  }

  if(!is.null(debiasing.weights) && length(debiasing.weights) != length(forest$Y.hat)) {
    stop("If specified, debiasing.weights must be a vector of length n.")
  }

  cluster.se <- length(forest$clusters) > 0

  clusters <- if (cluster.se) {
    forest$clusters
  } else {
    1:length(forest$Y.orig)
  }
  observation.weight <- observation_weights(forest)

  # Only use data selected via subsetting.
  subset.X.orig <- forest$X.orig[subset, , drop = FALSE]
  subset.W.orig <- forest$W.orig[subset]
  subset.W.hat <- forest$W.hat[subset]
  subset.Y.orig <- forest$Y.orig[subset]
  subset.Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]
  subset.clusters <- clusters[subset]
  subset.weights.raw <- observation.weight[subset]
  subset.weights <- subset.weights.raw / mean(subset.weights.raw)

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  binary.W <- all(subset.W.orig %in% c(0, 1))

  if (binary.W) {
    if (min(subset.W.hat) <= 0.01 && max(subset.W.hat) >= 0.99) {
      rng <- range(subset.W.hat)
      warning(paste0(
        "Estimated treatment propensities take values between ",
        round(rng[1], 3), " and ", round(rng[2], 3),
        " and in particular get very close to 0 and 1."
      ))
    }
  }

  if (!is.null(debiasing.weights)) {
    weights <- debiasing.weights[subset]
  } else if (binary.W) {
    IPW <- (subset.W.orig - subset.W.hat) / (subset.W.hat * (1 - subset.W.hat))
    weights <- IPW
  } else {
    # If treatment is continuous, estimate the variance of W given X
    # Details: https://github.com/grf-labs/grf/issues/608#issuecomment-582747602
    variance.forest <- regression_forest(subset.X.orig,
                                        (subset.W.orig - subset.W.hat)^2,
                                        clusters = subset.clusters,
                                        num.trees = num.trees.for.variance,
                                        ci.group.size = 1)
    V.hat <- predict(variance.forest)$predictions
    debiasing.weights <- subset.weights * (subset.W.orig - subset.W.hat) / V.hat
    weights <- debiasing.weights
  }

  # Compute doubly robust scores
  mu.w.hat <- subset.Y.hat + (subset.W.orig - subset.W.hat) * tau.hat.pointwise
  mu.w.hat.residual <- subset.Y.orig - mu.w.hat
  correction <- weights * mu.w.hat.residual
  Gamma.hat <- tau.hat.pointwise + correction

  if (!is.null(A)) {
    A <- as.matrix(A)
    if (nrow(A) == length(forest$Y.orig)) {
      A.subset <- A[subset, , drop = FALSE]
    } else if (nrow(A) == length(subset)) {
      A.subset <- A
    } else {
      stop("The number of rows of A does not match the number of training examples.")
    }
    if (is.null(colnames(A.subset))) {
      colnames(A.subset) <- paste0("A", 1:ncol(A))
    }
    DF <- data.frame(target = Gamma.hat, A.subset)
  } else {
    DF <- data.frame(target = Gamma.hat)
  }

  blp.ols <- lm(target ~ ., weights = subset.weights, data = DF)
  blp.summary <- lmtest::coeftest(blp.ols,
                                  vcov = sandwich::vcovCL,
                                  type = "HC3",
                                  cluster = subset.clusters
  )
  attr(blp.summary, "method") <-
    paste("Best linear projection of the conditional average treatment effect.",
          "Confidence intervals are cluster- and heteroskedasticity-robust (HC3)",
          sep="\n")

  blp.summary
}
