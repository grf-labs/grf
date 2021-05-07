#' Causal survival forest (experimental)
#'
#' Trains a causal survival forest that can be used to estimate
#' conditional average treatment effects tau(X). When the treatment assignment
#' is unconfounded, we have tau(X) = E[Y(1) - Y(0) | X = x],
#' where Y is the survival time up to a fixed maximum follow-up time.
#' Y(1) and Y(0) are potental outcomes corresponding to the two possible
#' treatment states.
#'
#'
#' An important assumption for identifying the conditional average treatment effect
#' tau(X) is that there exists a fixed positive constant M such that the probability
#' of observing an event time past the maximum follow-up time Y.max is at least M.
#' This may be an issue with data where most endpoint observations are censored.
#' The suggested resolution is to re-define the estimand as the treatment effect up
#' to some suitable maximum follow-up time Y.max. One can do this in practice by
#' thresholding Y before running causal_survival_forest: `D[Y >= Y.max] <- 1` and
#' `Y[Y >= Y.max] <- Y.max`. For details see Cui et al. (2020). The computational
#' complexity of this estimator scales with the cardinality of the event times Y.
#' If the number of samples is large and the Y grid dense, consider rounding the
#' event times (or supply a coarser grid with the `failure.times` argument).
#'
#' @param X The covariates.
#' @param Y The event time (may be negative).
#' @param W The treatment assignment (must be a binary vector with no NAs).
#' @param D The event type (0: censored, 1: failure).
#' @param W.hat Estimates of the treatment propensities E[W | Xi]. If W.hat = NULL,
#'              these are estimated using a separate regression forest. Default is NULL.
#' @param E1.hat Estimates of the expected survival time conditional on being treated
#'  E[Y | X = x, W = 1]. If E1.hat is NULL, then this is estimated with an S-learner using
#'  a survival forest.
#' @param E0.hat Estimates of the expected survival time conditional on being a control unit
#'  E[Y | X = x, W = 0]. If E0.hat is NULL, then this is estimated with an S-learner using
#'  a survival forest.
#' @param S.hat Estimates of the conditional survival function S(t, x, w) = P[Y > t | X = x, W = w].
#'  If S.hat is NULL, this is estimated using a survival forest. If provided:
#'  a N*T matrix of survival estimates. The grid should correspond to the T unique events in Y.
#'  Default is NULL.
#' @param C.hat Estimates of the conditional survival function for the censoring process S_C(t, x, w).
#'  If C.hat is NULL, this is estimated using a survival forest. If provided:
#'  a N*T matrix of survival estimates. The grid should correspond to the T unique events in Y.
#'  Default is NULL.
#' @param lambda.C.hat Estimates of the conditional hazard function -d/dt log(S_C(t, x, w)) for the censoring process.
#'  If lambda.C.hat is NULL, this is estimated from C.hat using a forward difference. If provided:
#'  a matrix of same dimensionality has C.hat.
#'  Default is NULL.
#' @param failure.times A vector of event times to fit the survival curves at. If NULL, then all the unique
#'  event times are used. This speeds up forest estimation by constraining the event grid. Observed event
#'  times are rounded down to the last sorted occurance less than or equal to the specified failure time.
#'  The time points should be in increasing order.
#'  Default is NULL.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param sample.weights Weights given to each sample in estimation.
#'                       If NULL, each observation receives the same weight.
#'                       Note: To avoid introducing confounding, weights should be
#'                       independent of the potential outcomes given X. Sample weights
#'                       are not used in survival spliting. Default is NULL.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'  Default is NULL (ignored).
#' @param equalize.cluster.weights If FALSE, each unit is given the same weight (so that bigger
#'  clusters get more weight). If TRUE, each cluster is given equal weight in the forest. In this case,
#'  during training, each tree uses the same number of observations from each drawn cluster: If the
#'  smallest cluster has K units, then when we sample a cluster during training, we only give a random
#'  K elements of the cluster to the tree-growing procedure. When estimating average treatment effects,
#'  each observation is given weight 1/cluster size, so that the total weight of each cluster is the
#'  same. Note that, if this argument is FALSE, sample weights may also be directly adjusted via the
#'  sample.weights argument. If this argument is TRUE, sample.weights must be set to NULL. Default is
#'  FALSE.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction. Default is 0.5.
#' @param mtry Number of variables tried for each split. Default is
#'             \eqn{\sqrt p + 20} where p is the number of variables.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#'                      Default is 5.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting). Default is TRUE.
#'  For a detailed description of honesty, honesty.fraction, honesty.prune.leaves, and recommendations for
#'  parameter tuning, see the grf algorithm reference.
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. Default is 0.5 (i.e. half of the data is used for
#'                         determining splits).
#' @param honesty.prune.leaves If TRUE, prunes the estimation sample tree such that no leaves
#'  are empty. If FALSE, keep the same tree as determined in the splits sample (if an empty leave is encountered, that
#'  tree is skipped and does not contribute to the estimate). Setting this to FALSE may improve performance on
#'  small/marginally powered data, but requires more trees (note: tuning does not adjust the number of trees).
#'  Only applies if honesty is enabled. Default is TRUE.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. Default is 0.05.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized. Default is 0.
#' @param stabilize.splits Whether or not the treatment and censoring status should be taken into account when
#'  determining the imbalance of a split. The requirement for valid split candidates is the same as in causal_forest
#'  with the additional constraint that num.failures(child) >= num.samples(parent) * alpha. Default is TRUE.
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2.
#' @param tune.parameters (Currently only applies to the regression forest used in W.hat estimation)
#'  A vector of parameter names to tune.
#'  If "all": all tunable parameters are tuned by cross-validation. The following parameters are
#'  tunable: ("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
#'   "honesty.prune.leaves", "alpha", "imbalance.penalty"). If honesty is FALSE the honesty.* parameters are not tuned.
#'  Default is "none" (no parameters are tuned).
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained causal_survival_forest forest object.
#'
#' @references Cui, Yifan, Michael R. Kosorok, Erik Sverdrup, Stefan Wager, and Ruoqing Zhu.
#'  "Estimating Heterogeneous Treatment Effects with Right-Censored Data via Causal Survival Forests."
#'  arXiv preprint arXiv:2001.09887, 2020.
#'
#' @examples
#' \donttest{
#' # Train a standard causal survival forest.
#' n <- 3000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y.max <- 1
#' failure.time <- pmin(rexp(n) * X[, 1] + W, Y.max)
#' censor.time <- 2 * runif(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' cs.forest <- causal_survival_forest(X, Y, W, D)
#'
#' # Predict using the forest.
#' X.test <- matrix(0.5, 10, p)
#' X.test[, 1] <- seq(0, 1, length.out = 10)
#' cs.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
#'
#' # Plot the estimated CATEs along with 95% confidence bands.
#' r.monte.carlo <- rexp(5000)
#' cate <- rep(NA, 10)
#' for (i in 1:10) {
#'   cate[i] <- mean(pmin(r.monte.carlo * X.test[i, 1] + 1, Y.max) -
#'                     pmin(r.monte.carlo * X.test[i, 1], Y.max))
#' }
#' plot(X.test[, 1], cate, type = 'l', col = 'red')
#' points(X.test[, 1], cs.pred$predictions)
#' lines(X.test[, 1], cs.pred$predictions + 2 * sqrt(cs.pred$variance.estimates), lty = 2)
#' lines(X.test[, 1], cs.pred$predictions - 2 * sqrt(cs.pred$variance.estimates), lty = 2)
#'
#' # Compute a doubly robust estimate of the average treatment effect.
#' average_treatment_effect(cs.forest)
#'
#' # Compute the best linear projection on the first covariate.
#' best_linear_projection(cs.forest, X[, 1])
#'
#' # Train the forest on a less granular grid.
#' cs.forest.grid <- causal_survival_forest(X, Y, W, D,
#'                                          failure.times = seq(min(Y), max(Y), length.out = 50))
#' plot(X.test[, 1], cs.pred$predictions)
#' points(X.test[, 1], predict(cs.forest.grid, X.test)$predictions, col = "blue")
#' }
#'
#' @export
causal_survival_forest <- function(X, Y, W, D,
                                   W.hat = NULL,
                                   E1.hat = NULL,
                                   E0.hat = NULL,
                                   S.hat = NULL,
                                   C.hat = NULL,
                                   lambda.C.hat = NULL,
                                   failure.times = NULL,
                                   num.trees = 2000,
                                   sample.weights = NULL,
                                   clusters = NULL,
                                   equalize.cluster.weights = FALSE,
                                   sample.fraction = 0.5,
                                   mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                                   min.node.size = 5,
                                   honesty = TRUE,
                                   honesty.fraction = 0.5,
                                   honesty.prune.leaves = TRUE,
                                   alpha = 0.05,
                                   imbalance.penalty = 0,
                                   stabilize.splits = TRUE,
                                   ci.group.size = 2,
                                   tune.parameters = "none",
                                   compute.oob.predictions = TRUE,
                                   num.threads = NULL,
                                   seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)
  D <- validate_observations(D, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)
  if (!all(W %in% c(0, 1))) {
    stop("The treatment values can only be 0 or 1.")
  }
  if (!all(D %in% c(0, 1))) {
    stop("The censor values can only be 0 or 1.")
  }
  if (sum(D) == 0) {
    stop("All observations are censored.")
  }
  if (is.null(failure.times)) {
    Y.grid <- sort(unique(Y))
  } else {
    if (is.unsorted(failure.times, strictly = TRUE)) {
      stop("Argument `failure.times` should be a vector with elements in increasing order.")
    }
    Y.grid <- failure.times
  }
  if (length(Y.grid) <= 2) {
    stop("The number of distinct event times should be more than 2.")
  }
  if (nrow(X) > 5000 && length(Y.grid) / nrow(X) > 0.1) {
    warning(paste0("The number of events are more than 10% of the sample size. ",
                   "To reduce the computational burden of fitting survival and ",
                   "censoring curves, consider rounding the event values `Y` or ",
                   "supplying a coarser grid with the `failure.times` argument. "))
  }

  args.orthog <- list(X = X,
                      Y = W,
                      num.trees = max(50, num.trees / 4),
                      sample.weights = sample.weights,
                      clusters = clusters,
                      equalize.cluster.weights = equalize.cluster.weights,
                      sample.fraction = sample.fraction,
                      mtry = mtry,
                      min.node.size = 5,
                      honesty = TRUE,
                      honesty.fraction = 0.5,
                      honesty.prune.leaves = TRUE,
                      alpha = alpha,
                      imbalance.penalty = imbalance.penalty,
                      ci.group.size = 1,
                      tune.parameters = tune.parameters,
                      compute.oob.predictions = TRUE,
                      num.threads = num.threads,
                      seed = seed)

  if (is.null(W.hat)) {
    forest.W <- do.call(regression_forest, args.orthog)
    W.hat <- predict(forest.W)$predictions
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, nrow(X))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }
  W.centered <- W - W.hat

  args.nuisance <- list(failure.times = failure.times,
                        num.trees = max(50, num.trees / 4),
                        sample.weights = sample.weights,
                        clusters = clusters,
                        equalize.cluster.weights = equalize.cluster.weights,
                        sample.fraction = sample.fraction,
                        mtry = mtry,
                        min.node.size = 15,
                        honesty = TRUE,
                        honesty.fraction = 0.5,
                        honesty.prune.leaves = TRUE,
                        alpha = alpha,
                        prediction.type = "Nelson-Aalen",
                        compute.oob.predictions = FALSE,
                        num.threads = num.threads,
                        seed = seed)

  # The survival function conditioning on being treated S(t, x, 1) estimated with an "S-learner".
  if (is.null(E1.hat)) {
    sf.survival <- do.call(survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.nuisance))
    S1.failure.times <- S0.failure.times <- sf.survival$failure.times
    # Computing OOB estimates for modified training samples is not a workflow we have implemented,
    # so we do it with a manual workaround here. Note that compute.oob.predictions has to be FALSE.
    sf.survival[["X.orig"]][, ncol(X) + 1] <- rep(1, nrow(X))
    S1.hat <- predict(sf.survival)$predictions
    sf.survival[["X.orig"]][, ncol(X) + 1] <- W
    E1.hat <- expected_survival(S1.hat, S1.failure.times)
  } else if (length(E1.hat) != nrow(X)) {
    stop("E1.hat has incorrect length.")
  }

  # The survival function conditioning on being a control unit S(t, x, 0) estimated with an "S-learner".
  if (is.null(E0.hat)) {
    if (!exists("sf.survival", inherits = FALSE)) {
      sf.survival <- do.call(survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.nuisance))
      S0.failure.times <- sf.survival$failure.times
    }
    sf.survival[["X.orig"]][, ncol(X) + 1] <- rep(0, nrow(X))
    S0.hat <- predict(sf.survival)$predictions
    sf.survival[["X.orig"]][, ncol(X) + 1] <- W
    E0.hat <- expected_survival(S0.hat, S0.failure.times)
  } else if (length(E0.hat) != nrow(X)) {
    stop("E0.hat has incorrect length.")
  }
  # Compute m(x) = e(X) E[T | X, W = 1] + (1 - e(X)) E[T | X, W = 0]
  m.hat <- W.hat * E1.hat + (1 - W.hat) * E0.hat

  # The conditional survival function S(t, x, w).
  if (is.null(S.hat)) {
    if (!exists("sf.survival", inherits = FALSE)) {
      sf.survival <- do.call(survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.nuisance))
    }
    # We want the predicted survival curves S.hat and C.hat on the common grid Y.grid:
    S.hat <- predict(sf.survival, failure.times = Y.grid)$predictions
  } else if (NROW(S.hat) != nrow(X)) {
    stop("S.hat has incorrect length.")
  } else if (NCOL(S.hat) != length(Y.grid)) {
    stop("S.hat has incorrect number of columns (should be equal to the number of events).")
  }

  # The conditional survival function for the censoring process S_C(t, x, w).
  if (is.null(C.hat)) {
    sf.censor <- do.call(survival_forest, c(list(X = cbind(X, W), Y = Y, D = 1 - D), args.nuisance))
    C.hat <- predict(sf.censor, failure.times = Y.grid)$predictions
  } else if (NROW(C.hat) != nrow(X)) {
    stop("C.hat has incorrect length.")
  } else if (NCOL(C.hat) != length(Y.grid)) {
    stop("C.hat has incorrect number of columns (should be equal to the number of events).")
  }

  if (any(C.hat == 0.0)) {
    stop("Some censoring probabilites are exactly zero.")
  }

  if (any(C.hat <= 0.05)) {
    warning(paste("Estimated censoring probabilites go as low as:", min(C.hat),
                "- an identifying assumption is that there exists a fixed positive constant M",
                "such that the probability of observing an event time past the maximum follow-up time Y.max",
                "is at least M. Formally, we assume: P(Y >= Y.max | X) > M.",
                "This warning appears when M is less than 0.05, at which point causal survival forest",
                "can not be expected to deliver reliable estimates."))
  } else if (any(C.hat < 0.2 & C.hat > 0.05)) {
    warning(paste("Estimated censoring probabilites are lower than 0.2.",
                  "An identifying assumption is that there exists a fixed positive constant M",
                  "such that the probability of observing an event time past the maximum follow-up time Y.max",
                  "is at least M. Formally, we assume: P(Y >= Y.max | X) > M."))
  }

  # Compute the pseudo outcomes
  eta <- compute_eta(S.hat, C.hat, lambda.C.hat, Y.grid, Y, D, m.hat, W.centered)
  validate_observations(eta[["numerator"]], X)
  validate_observations(eta[["denominator"]], X)

  data <- create_train_matrices(X,
                                treatment = W.centered,
                                survival.numerator = eta[["numerator"]],
                                survival.denominator = eta[["denominator"]],
                                censor = D,
                                sample.weights = sample.weights)

  args <- list(num.trees = num.trees,
               clusters = clusters,
               samples.per.cluster = samples.per.cluster,
               sample.fraction = sample.fraction,
               mtry = mtry,
               min.node.size = min.node.size,
               honesty = honesty,
               honesty.fraction = honesty.fraction,
               honesty.prune.leaves = honesty.prune.leaves,
               alpha = alpha,
               imbalance.penalty = imbalance.penalty,
               stabilize.splits = stabilize.splits,
               ci.group.size = ci.group.size,
               compute.oob.predictions = compute.oob.predictions,
               num.threads = num.threads,
               seed = seed)

  forest <- do.call.rcpp(causal_survival_train, c(data, args))
  class(forest) <- c("causal_survival_forest", "grf")
  forest[["eta"]] <- eta
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["D.orig"]] <- D
  forest[["W.hat"]] <- W.hat
  forest[["sample.weights"]] <- sample.weights
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["has.missing.values"]] <- has.missing.values

  forest
}

#' Predict with a causal survival forest forest
#'
#' Gets estimates of tau(X) using a trained causal survival forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions.
#'
#' @examples
#' \donttest{
#' # Train a standard causal survival forest.
#' n <- 3000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y.max <- 1
#' failure.time <- pmin(rexp(n) * X[, 1] + W, Y.max)
#' censor.time <- 2 * runif(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' cs.forest <- causal_survival_forest(X, Y, W, D)
#'
#' # Predict using the forest.
#' X.test <- matrix(0.5, 10, p)
#' X.test[, 1] <- seq(0, 1, length.out = 10)
#' cs.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
#'
#' # Plot the estimated CATEs along with 95% confidence bands.
#' r.monte.carlo <- rexp(5000)
#' cate <- rep(NA, 10)
#' for (i in 1:10) {
#'   cate[i] <- mean(pmin(r.monte.carlo * X.test[i, 1] + 1, Y.max) -
#'                     pmin(r.monte.carlo * X.test[i, 1], Y.max))
#' }
#' plot(X.test[, 1], cate, type = 'l', col = 'red')
#' points(X.test[, 1], cs.pred$predictions)
#' lines(X.test[, 1], cs.pred$predictions + 2 * sqrt(cs.pred$variance.estimates), lty = 2)
#' lines(X.test[, 1], cs.pred$predictions - 2 * sqrt(cs.pred$variance.estimates), lty = 2)
#'
#' # Compute a doubly robust estimate of the average treatment effect.
#' average_treatment_effect(cs.forest)
#'
#' # Compute the best linear projection on the first covariate.
#' best_linear_projection(cs.forest, X[, 1])
#'
#' # Train the forest on a less granular grid.
#' cs.forest.grid <- causal_survival_forest(X, Y, W, D,
#'                                          failure.times = seq(min(Y), max(Y), length.out = 50))
#' plot(X.test[, 1], cs.pred$predictions)
#' points(X.test[, 1], predict(cs.forest.grid, X.test)$predictions, col = "blue")
#' }
#'
#' @method predict causal_survival_forest
#' @export
predict.causal_survival_forest <- function(object,
                                           newdata = NULL,
                                           num.threads = NULL,
                                           estimate.variance = FALSE,
                                           ...) {
  # If possible, use pre-computed predictions.
  if (is.null(newdata) && !estimate.variance && !is.null(object$predictions)) {
    return(data.frame(predictions = object$predictions))
  }

  num.threads <- validate_num_threads(num.threads)

  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X)

  args <- list(forest.object = forest.short,
               num.threads = num.threads,
               estimate.variance = estimate.variance)

  if (!is.null(newdata)) {
    validate_newdata(newdata, X, allow.na = TRUE)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(causal_survival_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(causal_survival_predict_oob, c(train.data, args))
  }

  # Convert list to data frame.
  ret <- ret[c(1, 2)] # the last two entries are unused error estimates
  empty <- sapply(ret, function(elem) length(elem) == 0)
  do.call(cbind.data.frame, ret[!empty])
}

#' Compute pseudo outcomes (based on the influence function of tau(X) for each unit)
#' used for CART splitting.
#'
#' We compute:
#' eta = [sum (W - W.hat)^2 (1 / C.hat - integral_{0}^{Y_i} (lambdaC / C.hat) dt)]^(-1) ("denominator")
#'  * sum (Gamma / C.hat - integral_{0}^{Y_i} (lambdaC / C.hat)(W - W.hat)(Q(t) - m.hat)dt) ("numerator")
#'
#' where
#' Gamma = (W - W.hat)(Y - m.hat) if D == 1 else (W - W.hat)(Q(Y) - m.hat)
#' Q(t) - E[T | X, W, Y >= t]
#' m(x) = e(X) E[T | X, W = 1] + (1 - e(X)) E[T | X, W = 0]
#' lambda.C = -d/dt log(C.hat(t, x, w)) is the conditional hazard function of the censoring process.
#'
#' Some useful properties:
#' The expected survival time E[T] is the integral of the survival function S(t).
#' The conditional expected survival time E[T | T >= y] is y + the integral of S(t + y) / S(y).
#'
#' @param S.hat Estimates of the conditional survival function S(t, x, w).
#' @param C.hat Estimates of the conditional survival function for the censoring process S_C(t, x, w).
#' @param lambda.C.hat Estimates of the conditional hazard function for the censoring process S_C(t, x, w).
#' @param Y.grid The time values corresponding to S.hat and C.hat.
#' @param Y The event times.
#' @param D The censoring indicator.
#' @param m.hat Estimates of m(X).
#' @param W.centered W - W.hat.
#' @return A list with the numerator and denominator of eta.
#' @keywords internal
compute_eta <- function(S.hat,
                        C.hat,
                        lambda.C.hat,
                        Y.grid,
                        Y,
                        D,
                        m.hat,
                        W.centered) {
  # The event time values relabeled to consecutive integers 0/1 to length(Y.grid).
  Y.relabeled <- findInterval(Y, Y.grid)
  # S.hat and C.hat will have the same dimensions,
  # they are survival estimates corresponding to the same grid Y.grid.
  num.samples <- nrow(S.hat)
  grid.length <- ncol(S.hat)
  Y.diff <- diff(c(0, Y.grid))
  if (Y.diff[1] == 0) {
    Y.diff[1] = 1
  }

  # The censoring probabilities for the observed events.
  # Depending on the optional grid `failure.times` the smallest `Y.relabeled` may be zero,
  # which is why we subset appropriately here (for events before the first failure, C.Y.hat is one)
  C.Y.hat <- rep(1, num.samples)
  C.Y.hat[Y.relabeled != 0] <- C.hat[cbind(1:num.samples, Y.relabeled)]
  Y.relabeled[Y.relabeled == 0] = 1

  if (is.null(lambda.C.hat)) {
    # The conditional hazard function lambda.C.hat = -d/dt log(C.hat(t, x, w))
    # This simple forward difference approximation works reasonably well.
    log.surv.C <- -log(cbind(1, C.hat))
    lambda.C.hat <- log.surv.C[, 2:(grid.length + 1)] - log.surv.C[, 1:grid.length]
    lambda.C.hat <- sweep(lambda.C.hat, 2, Y.diff, "/")
  } else if (NROW(lambda.C.hat) != NROW(C.hat) || NCOL(lambda.C.hat) != NCOL(C.hat)) {
    stop("lambda.C.hat has incorrect dimensions.")
  }

  # The denominator simplifies to this
  denominator <- D * W.centered^2 / C.Y.hat

  # Even though we only need to compute the integral up to Yi, it is more clear (and as fast)
  # to simply compute everything for all time points t then at the end sum
  # up to the appropriate column Yi.

  # Compute Q(t, X) = E[T | X, W, T >= t]
  # We can quickly compute all these t conditional expectations by updating backwards.
  # For each time point t, the conditional expectation for sample i takes the form:
  # t + Y.diff[(t + 1):grid.length] %*% S.hat[i, t:(grid.length - 1)] / S.hat[i, t]
  Q.t.hat <- matrix(0, num.samples, grid.length)
  dot.products <- sweep(S.hat[, 1:(grid.length - 1)], 2, Y.diff[2:grid.length], "*")
  Q.t.hat[, 1] <- rowSums(dot.products)
  for (i in 2:(grid.length - 1)) {
    Q.t.hat[, i] <- Q.t.hat[, i - 1] - dot.products[, i - 1]
  }
  Q.t.hat <- Q.t.hat / S.hat
  Q.t.hat[is.infinite(Q.t.hat)] <- 0 # The points where S.hat = 0
  Q.t.hat <- sweep(Q.t.hat, 2, Y.grid, "+") # Add back t
  Q.t.hat[, grid.length] <- max(Y.grid)

  # Get Q.Y.hat = Q(Y, X) = E[T | X, W, T >= Y]
  Q.Y.hat <- Q.t.hat[cbind(1:num.samples, Y.relabeled)]
  numerator.one <- (D * (Y - m.hat) + (1 - D) * (Q.Y.hat - m.hat)) * W.centered / C.Y.hat

  integrand <- sweep(lambda.C.hat / C.hat * (Q.t.hat - m.hat), 2, Y.diff, "*")
  numerator.two <- rep(0, num.samples)
  # Store additional coefficients for computing robust scores later on
  integrand.update <- sweep(lambda.C.hat / C.hat, 2, Y.diff, "*")
  integral.update <- rep(0, num.samples)
  for (sample in 1:num.samples) {
    Y.index <- Y.relabeled[sample]
    numerator.two[sample] <- sum(integrand[sample, 1:Y.index]) * (W.centered[sample])
    integral.update[sample]<- sum(integrand.update[sample, 1:Y.index]) * (W.centered[sample])
  }

  numerator <- numerator.one - numerator.two

  list(numerator = numerator, denominator = denominator,
       numerator.one = numerator.one, numerator.two = numerator.two,
       integral.update = integral.update, C.Y.hat = C.Y.hat)
}

#' Compute E[T | X]
#'
#' @param S.hat The estimated survival curve.
#' @param Y.grid The time values corresponding to S.hat.
#' @return A vector of expected values.
#' @keywords internal
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))

  c(cbind(1, S.hat) %*% grid.diff)
}
