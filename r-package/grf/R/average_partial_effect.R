#' Estimate average partial effects using a causal forest
#' 
#' Gets estimates of the average partial effect,
#' 
#' The (conditional) average treatment effect (target.sample = all):
#'   1/n sum_{i = 1}^n Cov[Wi, Yi | X = Xi] / Var[Wi | X = Xi].
#'
#' Note that for a binary unconfounded treatment, the
#' average partial effect matches the average treatment effect.
#'
#' @param forest The trained forest.
#'
#' @examples \dontrun{
#' n = 2000; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 1/(1 + exp(-X[,2]))) + rnorm(n)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' tau.forest = causal_forest(X, Y, W)
#' tau.hat = predict(tau.forest, X.test)
#' average_partial_effect(tau.forest)
#' }
#'
#' @return An estimate of the average partial effect, along with standard error.
#' @export
average_partial_effect = function(forest) {

  if (!("causal_forest" %in% class(forest))) {
    stop("Average effect estimation only implemented for causal_forest")
  }
  
  if (is.null(forest$Y.hat) | is.null(forest$W.hat)) {
    stop("For average effect estimation to work, please train with precompute.nuisance = TRUE")
  }
  
  # This is a simple plugin estimate of the APE.
  tau.hat = predict(forest)$predictions
  cape.plugin = mean(tau.hat)
  
  # Estimate the variance of W given X. For binary treatments,
  # we get a good implicit estimator V.hat = e.hat (1 - e.hat), and
  # so this step is not needed. Note that if we use the present CAPE estimator
  # with a binary treatment and set V.hat = e.hat (1 - e.hat), then we recover
  # exactly the AIPW estimator of the CATE.
  variance_forest <- regression_forest(forest$X.orig,
                                       (forest$W.orig - forest$W.hat)^2)
  V.hat <- predict(variance_forest)$predictions
  
  # Compute a residual-based correction to the plugin.
  plugin.prediction = forest$Y.hat + (forest$W.orig - forest$W.hat) * tau.hat
  plugin.residual = forest$Y.orig - plugin.prediction
  cape.correction = mean((forest$W.orig - forest$W.hat) * plugin.residual / V.hat)
  cape.estimate = cape.plugin + cape.correction
  
  # Estimate variance using the calibration
  cape.se = sqrt(mean(((forest$W.orig - forest$W.hat) * plugin.residual / V.hat)^2) /
                   length(forest$W.orig))
    
  return(c(estimate=cape.estimate, std.err=cape.se))
}
