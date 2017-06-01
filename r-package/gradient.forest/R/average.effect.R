#' Estimate average treatment effects using a causal forest
#' 
#' Gets estimates of one of the following
#' 
#' The (conditional) average treatment effect (target.sample = all):
#'   sum_{i = 1}^n E[Y(1) - Y(0) | X = Xi] / n
#' The (conditional) average treatment effect on the treated (target.sample = treated):
#'   sum_{Wi = 1} E[Y(1) - Y(0) | X = Xi] / |{i : Wi = 1}|
#' The (conditional) average treatment effect on the controls (target.sample = control):
#'   sum_{Wi = 0} E[Y(1) - Y(0) | X = Xi] / |{i : Wi = 0}|
#'
#' @param forest The trained forest.
#' @param target.sample Which sample to aggregate treatment effects over.
#' @param method Method used for doubly robust inference. Can be either
#'               augmented inverse-propensity weighting (AIPW), or
#'               targeted maximum likelihood estimation (TMLE).
#'
#' @return Vector of predictions, along with (optional) variance estimates.
#' @export
estimate.average.effect = function(forest,
                                   target.sample=c("all", "treated", "control"),
                                   method=c("AIPW", "TMLE")) {

  if (class(forest) != "causal.forest") {
    stop("Average effect estimation only implemented for causal.forest")
  }
  
  if (is.null(forest$Y.hat) | is.null(forest$W.hat)) {
    stop("For average effect estimation to work, please train with precompute.nuisance = TRUE")
  }
  
  if (!all(forest$W.orig %in% c(0, 1))) {
    stop("Average effect estimation only implemented for binary treatment")
  }
  
  target.sample <- match.arg(target.sample)
  method <- match.arg(method)
  
  control.idx <- which(forest$W.orig == 0)
  treated.idx <- which(forest$W.orig == 1)
  
  # Retreive pointwise treatment effect predictions from forest, and
  # compute naive average effect estimates (notice that this uses OOB)
  tau.hat.pointwise <- predict(forest)$predictions
  if (target.sample == "all") {
    tau.avg.raw <- mean(tau.hat.pointwise)
  } else if (target.sample == "treated") {
    tau.avg.raw <- mean(tau.hat.pointwise[treated.idx])
  } else if (target.sample == "control") {
    tau.avg.raw <- mean(tau.hat.pointwise[control.idx])
  } else {
    stop("Invalid target sample.")
  }
  
  # Get estimates for the regress surfaces E[Y|X, W=0/1]
  Y.hat.0 <- forest$Y.hat - forest$W.hat * tau.hat.pointwise
  Y.hat.1 <- forest$Y.hat + (1 - forest$W.hat) * tau.hat.pointwise
  
  # Compute inverse-propensity-type weights that can be used for either
  # AIPW or TMLE estimation.
  
  if (target.sample == "all") {
    gamma.control.raw <- 1 / (1 - forest$W.hat[control.idx])
    gamma.treated.raw <- 1 / forest$W.hat[treated.idx]
  } else if (target.sample == "treated") {
    gamma.control.raw <- forest$W.hat[control.idx] / (1 - forest$W.hat[control.idx])
    gamma.treated.raw <- rep(1, length(treated.idx))
  } else if (target.sample == "control") {
    gamma.control.raw <- rep(1, length(control.idx))
    gamma.treated.raw <- (1 - forest$W.hat[treated.idx]) / forest$W.hat[treated.idx]
  } else {
    stop("Invalid target sample.")
  }
  
  gamma <- rep(0, length(forest$W.orig))
  gamma[control.idx] <- gamma.control.raw / sum(gamma.control.raw) * length(forest$W.orig)
  gamma[treated.idx] <- gamma.treated.raw / sum(gamma.treated.raw) * length(forest$W.orig)
  
  # Now compute a doubly robust correction
  if (method == "AIPW") {
    dr.correction.all <- forest$W.orig * gamma * (forest$Y.orig - Y.hat.1) -
      (1 - forest$W.orig) * gamma * (forest$Y.orig - Y.hat.0)
    dr.correction <- mean(dr.correction.all)
    sigma2.hat <- mean(dr.correction.all^2) / length(dr.correction.all)
  } else if (method == "TMLE") {
    stop("... not implemented ...")
  } else {
    stop("Invalid method.")
  }
  
  tau.avg <- tau.avg.raw + dr.correction
  tau.se <- sqrt(sigma2.hat)
  return(c(estimate=tau.avg, std.err=tau.se))
}