#' Estimate the average (conditional) local average treatment effect using a causal forest.
#'
#' Given an outcome Y, treatment W and instrument Z, the (conditional) local
#' average treatment effect is tau(x) = Cov[Y, Z | X = x] / Cov[W, Z | X = x].
#' This is the quantity that is estimated with an instrumental forest.
#' It can be intepreted causally in various ways. Given a homogeneity
#' assumption, tau(x) is simply the CATE at x. When W is binary
#' and there are no "defiers", Imbens and Angrist (1994) show that tau(x) can
#' be interpreted as an average treatment effect on compliers. This function
#' is about estimating tau = E[tau(X)] which, extending standard nomenclature,
#' should perhaps be called the Average (Conditional) Local Averate Treatment
#' Effect (ACLATE).
#' 
#' We estimate the ACLATE using a doubly robust estimator. See Chernozhukov
#' et al. (2016) for a discussion, and Section 5.2 of Athey and Wager (2017)
#' for an example using forests.
#'
#' If clusters are specified for the forest, then each cluster gets equal weight.
#' For example, if there are 10 clusters with 1 unit each and per-cluster ATE = 1,
#' and there are 10 clusters with 19 units each and per-cluster ATE = 0, then the
#' overall ATE is 0.5 (not 0.05).
#'
#' @param forest The trained forest.
#' @param compliance.score An estimate of the causal effect of Z on W,
#'                         i.e., Delta(X) = E[W | X, Z = 1] - E[W | X, Z = 0],
#'                         for each sample i = 1, ..., n.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the instrument Zi, treatment Wi or outcome Yi.
#'
#' @references Aronow, Peter M., and Allison Carnegie. "Beyond LATE: Estimation
#' of the average treatment effect with an instrumental variable." Political
#' Analysis 21.4 (2013): 492-506.
#' @references Athey, Susan, and Stefan Wager. "Efficient policy learning."
#' arXiv preprint arXiv:1702.02896 (2017).
#' @references Chernozhukov, Victor, Juan Carlos Escanciano, Hidehiko Ichimura,
#' Whitney K. Newey, and James M. Robins. "Locally robust semiparametric
#' estimation." arXiv preprint arXiv:1608.00033 (2016).
#' @references Imbens, Guido W., and Joshua D. Angrist. "Identification and
#' Estimation of Local Average Treatment Effects." Econometrica 62.2 (1994): 467-475.
#'
#' @return An estimate of the average (C)LATE, along with standard error.
#'
#' @importFrom stats coef lm predict var weighted.mean
#' @export
average_late <- function(forest,
                         compliance.score = NULL,
                         subset = NULL) {
  cluster.se <- length(forest$clusters) > 0

  if (!("instrumental_forest" %in% class(forest))) {
      stop(paste(
          "Average conditional local average treatment effect estimation",
          "only implemented for instrumental_forest."
      ))
  }
  
  if (!all(forest$Z.orig %in% c(0, 1))) {
      stop(paste(
          "Average conditional local average treatment effect estimation",
          "only implemented for binary instruments."
      ))
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

  clusters <- if (cluster.se) {
    forest$clusters
  } else {
    1:length(forest$Y.orig)
  }
  observation.weight <- observation_weights(forest)
  
  # The compliance forest estimates the effect of the "treatment" Z
  # on the "outcome" W.
  if (is.null(compliance.score)) {
      compliance.forest <- causal_forest(forest$X.orig,
                                         Y=forest$W.orig,
                                         W=forest$Z.orig,
                                         Y.hat=forest$W.hat,
                                         W.hat=forest$Z.hat,
                                         sample.weights = forest$sample.weights,
                                         clusters = clusters)
      compliance.score <- predict(compliance.forest)$predictions
  }

  # Only use data selected via subsetting.
  subset.W.orig <- forest$W.orig[subset]
  subset.W.hat <- forest$W.hat[subset]
  subset.Z.orig <- forest$Z.orig[subset]
  subset.Z.hat <- forest$Z.hat[subset]
  subset.Y.orig <- forest$Y.orig[subset]
  subset.Y.hat <- forest$Y.hat[subset]
  subset.compliance.score <- compliance.score[subset]
  subset.tau.hat.pointwise <- predict(forest)$predictions[subset]
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (min(subset.Z.hat) <= 0.01 || max(subset.Z.hat) >= 0.99) {
    rng <- range(subset.W.hat)
    warning(paste0(
      "Estimated treatment propensities take values between ",
      round(rng[1], 3), " and ", round(rng[2], 3),
      " and in particular get very close to 0 or 1. ",
      "Poor overlap may hurt perfmance for average conditional local average ",
      "treatment effect estimation."
    ))
  }
  
  if (min(subset.compliance.score) <= 0.01 * sd(subset.W.orig)) {
      warning(paste0(
          "The instrument appears to be weak, with some compliance scores as ",
          "low as ", round(min(subset.compliance.score), 4)
      ))
  }
  
  subset.g.hat = (subset.Z.orig - subset.Z.hat) /
      (subset.Z.hat * (1 - subset.Z.hat)) /
      subset.compliance.score
  Gamma.hat = subset.tau.hat.pointwise + subset.g.hat * 
      (subset.Y.orig - subset.Y.hat - 
       (subset.W.orig - subset.W.hat) * subset.tau.hat.pointwise)

  aclate.hat <- weighted.mean(Gamma.hat, subset.weights)
  
  if (cluster.se) {
      clust.avg <- Matrix::sparse.model.matrix(
          ~ factor(subset.clusters) + 0,
          transpose = TRUE
      ) %*% (Gamma.hat * subset.weights)
      sigma2.hat <- var(as.numeric(clust.avg)) * length(clust.avg) /
          sum(subset.weights)^2 
  } else {
      sigma2.hat <- var(Gamma.hat) / length(Gamma.hat)
  }
  
  return(c(estimate = aclate.hat, std.err = sqrt(sigma2.hat)))
}
