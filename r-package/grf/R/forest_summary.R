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
#' \donttest{
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
#' using a causal forest, or causal survival forest.
#'
#' Let tau(Xi) = E[Y(1) - Y(0) | X = Xi] be the CATE, and Ai be a vector of user-provided
#' covariates. This function provides a (doubly robust) fit to the linear model tau(Xi) ~ beta_0 + Ai * beta.
#'
#' Procedurally, we do so by regressing doubly robust scores derived from the
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
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained via the appropriate doubly robust score
#'               construction, e.g., in the case of causal_forests with a binary treatment, they
#'               are obtained via inverse-propensity weighting.
#' @param num.trees.for.weights In some cases (e.g., with causal forests with a continuous
#'               treatment), we need to train auxiliary forests to learn debiasing weights.
#'               This is the number of trees used for this task. Note: this argument is only
#'               used when debiasing.weights = NULL.
#'
#' @references Chernozhukov, Victor, and Vira Semenova. "Simultaneous inference for
#'             Best Linear Predictor of the Conditional Average Treatment Effect and
#'             other structural functions." arXiv preprint arXiv:1702.06240 (2017).
#'
#' @examples
#' \donttest{
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
#' @export
best_linear_projection <- function(forest,
                                   A = NULL,
                                   subset = NULL,
                                   debiasing.weights = NULL,
                                   num.trees.for.weights = 500) {
  clusters <- if (length(forest$clusters) > 0) {
    forest$clusters
  } else {
    1:length(forest$Y.orig)
  }
  observation.weight <- observation_weights(forest)

  subset <- validate_subset(forest, subset)
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  if (!is.null(debiasing.weights)) {
    if (length(debiasing.weights) == length(forest$Y.orig)) {
      debiasing.weights <- debiasing.weights[subset]
    } else if (length(debiasing.weights) != length(subset)) {
      stop("If specified, debiasing.weights must be a vector of length n or the subset length.")
    }
  }

  binary.W <- all(forest$W.orig %in% c(0, 1))

  if (binary.W) {
    if (min(forest$W.hat[subset]) <= 0.01 && max(forest$W.hat[subset]) >= 0.99) {
      rng <- range(forest$W.hat[subset])
      warning(paste0(
        "Estimated treatment propensities take values between ",
        round(rng[1], 3), " and ", round(rng[2], 3),
        " and in particular get very close to 0 and 1."
      ))
    }
  }

  if (any(c("causal_forest", "causal_survival_forest") %in% class(forest))) {
    DR.scores <- get_scores(forest, subset = subset, debiasing.weights = debiasing.weights,
                            num.trees.for.weights = num.trees.for.weights)
  } else {
    stop("`best_linear_projection` is only implemented for `causal_forest` and `causal_survival_forest`")
  }

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
    DF <- data.frame(target = DR.scores, A.subset)
  } else {
    DF <- data.frame(target = DR.scores)
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
