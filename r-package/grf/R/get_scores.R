#' Compute doubly robust (AIPW) scores for average treatment effect estimation
#' using a causal forest. Under regularity conditions, the average of the DR.scores
#' is an efficient estimate of the average treatment effect.
#'
#' @param forest A trained causal forest.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) they are obtained via inverse-propensity weighting.
#' 
#' @export
get_scores_ATE = function(forest,
                          subset,
                          debiasing.weights = NULL) {

  if (!("causal_forest" %in% class(forest))) {
    stop("The forest must be a causal_forest.")
  }

  if(!all(forest$W.orig %in% c(0, 1))) {
    stop("ATE is only implemented for a binary treatment. See APE for continuous case.")
  }

  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]

  # Form AIPW scores. Note: We are implicitly using the following
  # estimates for the regression surfaces E[Y|X, W=0/1]:
  # Y.hat.0 <- Y.hat - W.hat * tau.hat.pointwise
  # Y.hat.1 <- Y.hat + (1 - W.hat) * tau.hat.pointwise

  Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))
  if (is.null(debiasing.weights)) {
    debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
  }
  tau.hat.pointwise + debiasing.weights * Y.residual
}

#' Compute doubly robust scores for average partial effect estimation with continuous
#' treatment, using a causal forest.
#'
#' @param forest A trained causal forest.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained by estimating Var[W | X = x] using a new forest.
#' @param num.trees.for.weights Number of trees used to estimate Var[W | X = x]. Note: this
#'               argument is only used when debiasing.weights = NULL.
#'               
#' @export
get_scores_APE = function(forest,
                          subset,
                          debiasing.weights = NULL,
                          num.trees.for.weights = 500) {
  if (!("causal_forest" %in% class(forest))) {
    stop("The forest must be a causal_forest.")
  }
  
  # Start by learning debiasing weights if needed.
  # The goal is to estimate the variance of W given X. For binary treatments,
  # we get a good implicit estimator V.hat = e.hat (1 - e.hat), and
  # so this step is not needed. Note that if we use the present CAPE estimator
  # with a binary treatment and set V.hat = e.hat (1 - e.hat), then we recover
  # exactly the AIPW estimator of the CATE.
  if (is.null(debiasing.weights)) {
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
                                         ci.group.size = 1)
    V.hat <- predict(variance_forest)$predictions
    debiasing.weights.all <- (forest$W.orig - forest$W.hat) / V.hat
    debiasing.weights <- debiasing.weights.all[subset]
  }
  
  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]
  
  # Form AIPW-type scores.
  Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))
  tau.hat.pointwise + debiasing.weights * Y.residual
}