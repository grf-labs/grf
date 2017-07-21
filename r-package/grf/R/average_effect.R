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
#' @examples
#' # Train a causal forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' c.forest = causal_forest(X, Y, W)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#' # Estimate the conditional average treatment effect on the full sample (CATE).
#' estimate_average_effect(c.forest, target.sample = "all")
#' 
#' # Estimate the conditional average treatment effect on the treated sample (CATT).
#' # We don't expect much difference between the CATE and the CATT in this example,
#' # since treatment assignment was randomized.
#' estimate_average_effect(c.forest, target.sample = "treated")
#'
#' @return Vector of predictions, along with (optional) variance estimates.
#' @export
estimate_average_effect = function(forest,
                                   target.sample=c("all", "treated", "control"),
                                   method=c("AIPW", "TMLE")) {

  if (!("causal_forest" %in% class(forest))) {
    stop("Average effect estimation only implemented for causal_forest")
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
  
  if (method == "TMLE") {
    loaded <- requireNamespace("sandwich", quietly = TRUE)
    if (!loaded) {
      warning("To use TMLE, please install the package `sandwich`. Using AIPW instead.")
      method = "AIPW"
    }
  }
  
  # Now apply a doubly robust correction
  if (method == "AIPW") {
  
    # Compute normalized inverse-propensity-type weights gamma
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
  
    dr.correction.all <- forest$W.orig * gamma * (forest$Y.orig - Y.hat.1) -
      (1 - forest$W.orig) * gamma * (forest$Y.orig - Y.hat.0)
    dr.correction <- mean(dr.correction.all)
    sigma2.hat <- mean(dr.correction.all^2) / length(dr.correction.all)
    
  } else if (method == "TMLE") {
    
    if (target.sample == "all") {
      eps.tmle.robust.0 <-
        lm(B ~ A + 0, data=data.frame(A=1/(1 - forest$W.hat[forest$W.orig==0]),
                                      B=forest$Y.orig[forest$W.orig==0]-Y.hat.0[forest$W.orig==0]))
      eps.tmle.robust.1 <-
        lm(B ~ A + 0, data=data.frame(A=1/forest$W.hat[forest$W.orig==1],
                                      B=forest$Y.orig[forest$W.orig==1]-Y.hat.1[forest$W.orig==1]))
      delta.tmle.robust.0 <- predict(eps.tmle.robust.0, newdata=data.frame(A=mean(1/(1 - forest$W.hat))))
      delta.tmle.robust.1 <- predict(eps.tmle.robust.1, newdata=data.frame(A=mean(1/forest$W.hat)))
      dr.correction <- delta.tmle.robust.1 - delta.tmle.robust.0
      # use robust SE
      sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * mean(1/(1 - forest$W.hat))^2 +
        sandwich::vcovHC(eps.tmle.robust.1) * mean(1/forest$W.hat)^2
    } else if (target.sample == "treated") {
      eps.tmle.robust.0 <-
        lm(B ~ A + 0,
           data=data.frame(A=forest$W.hat[forest$W.orig==0]/(1 - forest$W.hat[forest$W.orig==0]),
                           B=forest$Y.orig[forest$W.orig==0]-Y.hat.0[forest$W.orig==0]))
      new.center <- mean(forest$W.hat[forest$W.orig==1]/(1 - forest$W.hat[forest$W.orig==1]))
      delta.tmle.robust.0 <- predict(eps.tmle.robust.0,
                                     newdata=data.frame(A=new.center))
      dr.correction <- -delta.tmle.robust.0
      sigma2.hat = sandwich::vcovHC(eps.tmle.robust.0) * new.center^2 +
        var(forest$Y.orig[forest$W.orig==1]-Y.hat.1[forest$W.orig==1]) / sum(forest$W.orig==1)
    } else if (target.sample == "control") {
      eps.tmle.robust.1 <-
        lm(B ~ A + 0,
           data=data.frame(A=(1 - forest$W.hat[forest$W.orig==1])/forest$W.hat[forest$W.orig==1],
                           B=forest$Y.orig[forest$W.orig==1]-Y.hat.1[forest$W.orig==1]))
      new.center <- mean((1 - forest$W.hat[forest$W.orig==0])/forest$W.hat[forest$W.orig==0])
      delta.tmle.robust.1 <- predict(eps.tmle.robust.1,
                                     newdata=data.frame(A=new.center))
      dr.correction <- delta.tmle.robust.1
      sigma2.hat = var(forest$Y.orig[forest$W.orig==0]-Y.hat.0[forest$W.orig==0]) / sum(forest$W.orig==0) +
        sandwich::vcovHC(eps.tmle.robust.1) * new.center^2
    } else {
      stop("Invalid target sample.")
    }
    
  } else {
    stop("Invalid method.")
  }
  
  tau.avg <- tau.avg.raw + dr.correction
  tau.se <- sqrt(sigma2.hat)
  return(c(estimate=tau.avg, std.err=tau.se))
}
