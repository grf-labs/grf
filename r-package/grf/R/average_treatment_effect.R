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
#' @param forest The trained forest.
#' @param target.sample Which sample to aggregate treatment effects over.
#' @param method Method used for doubly robust inference. Can be either
#'               augmented inverse-propensity weighting (AIPW), or
#'               targeted maximum likelihood estimation (TMLE).
#'
#' @examples \dontrun{
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
#' average_treatment_effect(c.forest, target.sample = "all")
#' 
#' # Estimate the conditional average treatment effect on the treated sample (CATT).
#' # We don't expect much difference between the CATE and the CATT in this example,
#' # since treatment assignment was randomized.
#' average_treatment_effect(c.forest, target.sample = "treated")
#' }
#'
#' @return An estimate of the average treatment effect, along with standard error.
#' @export
average_treatment_effect = function(forest,
                                    target.sample=c("all", "treated", "control", "overlap"),
                                    method=c("AIPW", "TMLE")) {

  target.sample <- match.arg(target.sample)
  method <- match.arg(method)
  cluster.se <- length(forest$clusters) > 0
  
  if (!("causal_forest" %in% class(forest))) {
    stop("Average effect estimation only implemented for causal_forest")
  }
  
  if (is.null(forest$Y.hat) | is.null(forest$W.hat)) {
    stop("For average effect estimation to work, please train with precompute.nuisance = TRUE")
  }
  
  #
  # Address the overlap case separately, as this is a very different estimation problem.
  # The method argument (AIPW vs TMLE) is ignored in this case, as both methods are effectively
  # the same here. Also, note that the overlap-weighted estimator generalizes naturally to the
  # non-binary W case -- see, e.g., Robinson (Econometrica, 1988) -- and so we do not require
  # W to be binary here.
  #
  
  if (target.sample == "overlap") {
    W.residual <- forest$W.orig - forest$W.hat
    Y.residual <- forest$Y.orig - forest$Y.hat
    tau.ols <- lm(Y.residual ~ W.residual)
    tau.est <- coef(summary(tau.ols))[2,1]
    
    if (cluster.se) {
      tau.se <- sqrt(sandwich::vcovCL(tau.ols, cluster = forest$clusters)[2,2])
    } else {
      tau.se <- sqrt(sandwich::vcovHC(tau.ols)[2,2])
    }
    
    return(c(estimate=tau.est, std.err=tau.se))
  }
  
  if (!all(forest$W.orig %in% c(0, 1))) {
    stop(paste("Average treatment effect estimation only implemented for binary treatment.",
               "See `average_partial_effect` for continuous W."))
  }
  
  if (min(forest$W.hat) <= 0.01 && max(forest$W.hat) >= 0.99) {
    rng = range(forest$W.hat)
    warning(paste0("Estimated treatment propensities take values between ",
                   round(rng[1], 3), " and ", round(rng[2], 3),
                   " and in particular get very close to 0 and 1. ",
                   "In this case, using `target.sample=overlap`, or filtering data as in ",
                   "Crump, Hotz, Imbens, and Mitnik (Biometrika, 2009) may be helpful."))
  } else if (min(forest$W.hat) <= 0.01 && target.sample != "treated") {
    warning(paste0("Estimated treatment propensities go as low as ",
                   round(min(forest$W.hat), 3), " which means that treatment ",
                   "effects for some controls may not be well identified. ",
                   "In this case, using `target.sample=treated` may be helpful."))
  } else if (max(forest$W.hat) >= 0.99 && target.sample != "control") {
    warning(paste0("Estimated treatment propensities go as high as ",
                   round(max(forest$W.hat), 3), " which means that treatment ",
                   "effects for some treated units may not be well identified. ",
                   "In this case, using `target.sample=control` may be helpful."))
  }
  
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
    
    if (cluster.se) {
      correction.clust <- Matrix::sparse.model.matrix(
        ~ factor(forest$clusters) + 0,
        transpose = TRUE) %*% dr.correction.all
      sigma2.hat <- sum(correction.clust^2) / length(dr.correction.all) /
        (length(dr.correction.all) - 1)
    } else {
      sigma2.hat <- mean(dr.correction.all^2) / (length(dr.correction.all) - 1)
    }
    
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
      if (cluster.se) {
        sigma2.hat <- sandwich::vcovCL(eps.tmle.robust.0, cluster = forest$clusters[forest$W.orig==0]) *
          mean(1/(1 - forest$W.hat))^2 +
          sandwich::vcovCL(eps.tmle.robust.1, cluster = forest$clusters[forest$W.orig==1]) *
          mean(1/forest$W.hat)^2
      } else {
        sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * mean(1/(1 - forest$W.hat))^2 +
          sandwich::vcovHC(eps.tmle.robust.1) * mean(1/forest$W.hat)^2
      }
    } else if (target.sample == "treated") {
      eps.tmle.robust.0 <-
        lm(B ~ A + 0,
           data=data.frame(A=forest$W.hat[forest$W.orig==0]/(1 - forest$W.hat[forest$W.orig==0]),
                           B=forest$Y.orig[forest$W.orig==0]-Y.hat.0[forest$W.orig==0]))
      new.center <- mean(forest$W.hat[forest$W.orig==1]/(1 - forest$W.hat[forest$W.orig==1]))
      delta.tmle.robust.0 <- predict(eps.tmle.robust.0,
                                     newdata=data.frame(A=new.center))
      dr.correction <- -delta.tmle.robust.0
      if (cluster.se) {
        s.0 <- sandwich::vcovCL(eps.tmle.robust.0, cluster = forest$clusters[forest$W.orig==0]) *
          new.center^2
        delta.1 <- Matrix::sparse.model.matrix(
          ~ factor(forest$clusters[forest$W.orig==1]) + 0,
          transpose = TRUE) %*% (forest$Y.orig[forest$W.orig==1]-Y.hat.1[forest$W.orig==1])
        s.1 <- sum(delta.1^2) / sum(forest$W.orig==1) / (sum(forest$W.orig==1) - 1)
        sigma2.hat <- s.0 + s.1
      } else {
        sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * new.center^2 +
          var(forest$Y.orig[forest$W.orig==1]-Y.hat.1[forest$W.orig==1]) / sum(forest$W.orig==1)
      }
    } else if (target.sample == "control") {
      eps.tmle.robust.1 <-
        lm(B ~ A + 0,
           data=data.frame(A=(1 - forest$W.hat[forest$W.orig==1])/forest$W.hat[forest$W.orig==1],
                           B=forest$Y.orig[forest$W.orig==1]-Y.hat.1[forest$W.orig==1]))
      new.center <- mean((1 - forest$W.hat[forest$W.orig==0])/forest$W.hat[forest$W.orig==0])
      delta.tmle.robust.1 <- predict(eps.tmle.robust.1,
                                     newdata=data.frame(A=new.center))
      dr.correction <- delta.tmle.robust.1
      if (cluster.se) {
        delta.0 <- Matrix::sparse.model.matrix(
          ~ factor(forest$clusters[forest$W.orig==0]) + 0,
          transpose = TRUE) %*% (forest$Y.orig[forest$W.orig==0]-Y.hat.0[forest$W.orig==0])
        s.0 <- sum(delta.0^2) / sum(forest$W.orig==0) / (sum(forest$W.orig==0) - 1)
        s.1 <- sandwich::vcovCL(eps.tmle.robust.1, cluster = forest$clusters[forest$W.orig==1]) *
          new.center^2
        sigma2.hat <- s.0 + s.1
      } else {
        sigma2.hat <- var(forest$Y.orig[forest$W.orig==0]-Y.hat.0[forest$W.orig==0]) / sum(forest$W.orig==0) +
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
  return(c(estimate=tau.avg, std.err=tau.se))
}
