#' Causal survival forest
#'
#' Trains a causal survival forest that can be used to estimate
#' conditional treatment effects tau(X) with right-censored outcomes.
#' We estimate either 1)
#' tau(X) = E[min(T(1), horizon) - min(T(0), horizon) | X = x],
#' where T(1) and T(0) are potental outcomes corresponding to the two possible treatment states
#' and `horizon` is the maximum follow-up time, or 2)
#' tau(X) = P(T(1) > horizon | X = x) - P(T(0) > horizon | X = x), for a chosen time point `horizon`.
#'
#'
#' The causal survival forest paper defines the survival function in the 2nd estimand with weak inequality.
#' It is defined using strict inequality in the R package (note that P(T >= h) = P(T > h - epsilon)).
#'
#' @param X The covariates.
#' @param Y The event time (must be non-negative).
#' @param W The treatment assignment (must be a binary vector with no NAs).
#' @param D The event type (0: censored, 1: failure).
#' @param W.hat Estimates of the treatment propensities E[W | X = x]. If W.hat = NULL,
#'              these are estimated using a separate regression forest. Default is NULL.
#' @param target The target estimand. Choices are Restricted Mean Survival Time ("RMST") which estimates 1)
#'  E[min(T(1), horizon) - min(T(0), horizon) | X = x], or "survival.probability" which estimates 2)
#'  P(T(1) > horizon | X = x) - P(T(0) > horizon | X = x). Default is "RMST".
#' @param horizon A scalar that defines the estimand (required). If target is "RMST" then this defines
#'  the maximum follow-up time. If target is "survival.probability", then this defines the time point
#'  for the absolute risk difference estimate.
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
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. This parameter plays the same
#'  role as in causal forest and survival forest, where for the latter the number of failures in
#'  each child has to be at least one or `alpha` times the number of samples in the parent node. Default is 0.05.
#'  (On data with very low event rate the default value may be too high for the forest to split
#'  and lowering it may be beneficial).
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
#'  Journal of the Royal Statistical Society: Series B, forthcoming.
#'
#' @examples
#' \donttest{
#' # Train a causal survival forest targeting a Restricted Mean Survival Time (RMST)
#' # with maximum follow-up time set to `horizon`.
#' n <- 2000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' horizon <- 1
#' failure.time <- pmin(rexp(n) * X[, 1] + W, horizon)
#' censor.time <- 2 * runif(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' # Save computation time by constraining the event grid by discretizing (rounding) continuous events.
#' cs.forest <- causal_survival_forest(X, round(Y, 2), W, D, horizon = horizon)
#' # Or do so more flexibly by defining your own time grid using the failure.times argument.
#' # grid <- seq(min(Y), max(Y), length.out = 150)
#' # cs.forest <- causal_survival_forest(X, Y, W, D, horizon = horizon, failure.times = grid)
#'
#' # Predict using the forest.
#' X.test <- matrix(0.5, 10, p)
#' X.test[, 1] <- seq(0, 1, length.out = 10)
#' cs.pred <- predict(cs.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' cs.pred <- predict(cs.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' c.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
#'
#' # Compute a doubly robust estimate of the average treatment effect.
#' average_treatment_effect(cs.forest)
#'
#' # Compute the best linear projection on the first covariate.
#' best_linear_projection(cs.forest, X[, 1])
#'
#' # See if a causal survival forest succeeded in capturing heterogeneity by plotting
#' # the TOC and calculating a 95% CI for the AUTOC.
#' train <- sample(1:n, n / 2)
#' eval <- -train
#' train.forest <- causal_survival_forest(X[train, ], Y[train], W[train], D[train], horizon = horizon)
#' eval.forest <- causal_survival_forest(X[eval, ], Y[eval], W[eval], D[eval], horizon = horizon)
#' rate <- rank_average_treatment_effect(eval.forest,
#'                                       predict(train.forest, X[eval, ])$predictions)
#' plot(rate)
#' paste("AUTOC:", round(rate$estimate, 2), "+/", round(1.96 * rate$std.err, 2))
#' }
#'
#' @export
causal_survival_forest <- function(X, Y, W, D,
                                   W.hat = NULL,
                                   target = c("RMST", "survival.probability"),
                                   horizon = NULL,
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
  target <- match.arg(target)
  if (is.null(horizon) || !is.numeric(horizon) || length(horizon) != 1) {
    stop("The `horizon` argument defining the estimand is required.")
  }
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)
  D <- validate_observations(D, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)
  if (any(Y < 0)) {
    stop("The event times must be non-negative.")
  }
  if (!all(W %in% c(0, 1))) {
    stop("The treatment values can only be 0 or 1.")
  }
  if (!all(D %in% c(0, 1))) {
    stop("The censor values can only be 0 or 1.")
  }
  if (sum(D) == 0) {
    stop("All observations are censored.")
  }
  if (target == "RMST") {
    # f(T) <- min(T, horizon)
    D[Y >= horizon] <- 1
    Y[Y >= horizon] <- horizon
    fY <- Y
  } else {
    # f(T) <- 1{T > horizon}
    fY <- as.numeric(Y > horizon)
  }
  if (is.null(failure.times)) {
    Y.grid <- sort(unique(Y))
  } else if (min(Y) < min(failure.times)) {
    stop("If provided, `failure.times` should be a grid starting on or before min(Y).")
  } else {
    Y.grid <- failure.times
  }
  if (length(Y.grid) <= 2) {
    stop("The number of distinct event times should be more than 2.")
  }
  if (horizon < min(Y.grid)) {
    stop("`horizon` cannot be before the first event.")
  }
  if (nrow(X) > 5000 && length(Y.grid) / nrow(X) > 0.1) {
    warning(paste0("The number of events are more than 10% of the sample size. ",
                   "To reduce the computational burden of fitting survival and ",
                   "censoring curves, consider discretizing the event values `Y` or ",
                   "supplying a coarser grid with the `failure.times` argument. "), immediate. = TRUE)
  }

  if (is.null(W.hat)) {
    forest.W <- regression_forest(X, W, num.trees = max(50, num.trees / 4),
                                  sample.weights = sample.weights, clusters = clusters,
                                  equalize.cluster.weights = equalize.cluster.weights,
                                  sample.fraction = sample.fraction, mtry = mtry,
                                  min.node.size = 5, honesty = TRUE,
                                  honesty.fraction = 0.5, honesty.prune.leaves = TRUE,
                                  alpha = alpha, imbalance.penalty = imbalance.penalty,
                                  ci.group.size = 1, tune.parameters = tune.parameters,
                                  compute.oob.predictions = TRUE,
                                  num.threads = num.threads, seed = seed)
    W.hat <- predict(forest.W)$predictions
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, nrow(X))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }
  W.centered <- W - W.hat

  args.nuisance <- list(failure.times = failure.times,
                        num.trees = max(50, min(num.trees / 4, 500)),
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
                        prediction.type = "Nelson-Aalen", # to guarantee non-zero estimates.
                        compute.oob.predictions = FALSE,
                        num.threads = num.threads,
                        seed = seed)

  # E[f(T) | X] = e(X) E[f(T) | X, W = 1] + (1 - e(X)) E[f(T) | X, W = 0]
  sf.survival <- do.call(survival_forest, c(list(X = cbind(X, W), Y = Y, D = D), args.nuisance))
  # The survival function conditioning on being treated S(t, x, 1) estimated with an "S-learner".
  # Computing OOB estimates for modified training samples is not a workflow we have implemented,
  # so we do it with a manual workaround here. Note that compute.oob.predictions has to be FALSE.
  sf.survival[["X.orig"]][, ncol(X) + 1] <- rep(1, nrow(X))
  S1.hat <- predict(sf.survival)$predictions
  # The survival function conditioning on being a control unit S(t, x, 0) estimated with an "S-learner".
  sf.survival[["X.orig"]][, ncol(X) + 1] <- rep(0, nrow(X))
  S0.hat <- predict(sf.survival)$predictions
  sf.survival[["X.orig"]][, ncol(X) + 1] <- W
  if (target == "RMST") {
    Y.hat <- W.hat * expected_survival(S1.hat, sf.survival$failure.times) +
      (1 - W.hat) * expected_survival(S0.hat, sf.survival$failure.times)
  } else {
    horizonS.index <- findInterval(horizon, sf.survival$failure.times)
    if (horizonS.index == 0) {
      Y.hat <- rep(1, nrow(X))
    } else {
      Y.hat <- W.hat * S1.hat[, horizonS.index] + (1 - W.hat) * S0.hat[, horizonS.index]
    }
  }

  # The conditional survival function S(t, x, w).
  S.hat <- predict(sf.survival, failure.times = Y.grid)$predictions
  # The conditional survival function for the censoring process S_C(t, x, w).
  args.nuisance$compute.oob.predictions <- TRUE
  sf.censor <- do.call(survival_forest, c(list(X = cbind(X, W), Y = Y, D = 1 - D), args.nuisance))
  C.hat <- predict(sf.censor, failure.times = Y.grid)$predictions
  if (target == "survival.probability") {
    # Evaluate psi up to horizon
    D[Y > horizon] <- 1
    Y[Y > horizon] <- horizon
  }

  Y.index <- findInterval(Y, Y.grid) # (invariance: Y.index > 0)
  C.Y.hat <- C.hat[cbind(seq_along(Y.index), Y.index)] # Pick out P[Ci > Yi | Xi, Wi]

  if (target == "RMST" && any(C.Y.hat <= 0.05)) {
    warning(paste("Estimated censoring probabilities go as low as:", round(min(C.Y.hat), 5),
                  "- an identifying assumption is that there exists a fixed positive constant M",
                  "such that the probability of observing an event past the maximum follow-up time ",
                  "is at least M (i.e. P(T > horizon | X) > M).",
                  "This warning appears when M is less than 0.05, at which point causal survival forest",
                  "can not be expected to deliver reliable estimates."), immediate. = TRUE)
  } else if (target == "RMST" && any(C.Y.hat < 0.2)) {
    warning(paste("Estimated censoring probabilities are lower than 0.2",
                  "- an identifying assumption is that there exists a fixed positive constant M",
                  "such that the probability of observing an event past the maximum follow-up time ",
                  "is at least M (i.e. P(T > horizon | X) > M)."))
  } else if (target == "survival.probability" && any(C.Y.hat <= 0.001)) {
    warning(paste("Estimated censoring probabilities go as low as:", round(min(C.Y.hat), 5),
                  "- forest estimates will likely be very unstable, a larger target `horizon`",
                  "is recommended."), immediate. = TRUE)
  } else if (target == "survival.probability" && any(C.Y.hat < 0.05)) {
    warning(paste("Estimated censoring probabilities are lower than 0.05",
                  "and forest estimates may not be stable. Using a smaller target `horizon`",
                  "may help."))
  }

  psi <- compute_psi(S.hat, C.hat, C.Y.hat, Y.hat, W.centered,
                     D, fY, Y.index, Y.grid, target, horizon)
  validate_observations(psi[["numerator"]], X)
  validate_observations(psi[["denominator"]], X)

  data <- create_train_matrices(X,
                                treatment = W.centered,
                                survival.numerator = psi[["numerator"]],
                                survival.denominator = psi[["denominator"]],
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
  forest[["seed"]] <- seed
  forest[["_psi"]] <- psi
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["D.orig"]] <- D
  forest[["Y.hat"]] <- Y.hat
  forest[["W.hat"]] <- W.hat
  forest[["sample.weights"]] <- sample.weights
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["has.missing.values"]] <- has.missing.values
  forest[["target"]] <- target
  forest[["horizon"]] <- horizon

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
#' @return Vector of predictions along with optional variance estimates.
#'
#' @examples
#' \donttest{
#' # Train a causal survival forest targeting a Restricted Mean Survival Time (RMST)
#' # with maximum follow-up time set to `horizon`.
#' n <- 2000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' horizon <- 1
#' failure.time <- pmin(rexp(n) * X[, 1] + W, horizon)
#' censor.time <- 2 * runif(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' # Save computation time by constraining the event grid by discretizing (rounding) continuous events.
#' cs.forest <- causal_survival_forest(X, round(Y, 2), W, D, horizon = horizon)
#' # Or do so more flexibly by defining your own time grid using the failure.times argument.
#' # grid <- seq(min(Y), max(Y), length.out = 150)
#' # cs.forest <- causal_survival_forest(X, Y, W, D, horizon = horizon, failure.times = grid)
#'
#' # Predict using the forest.
#' X.test <- matrix(0.5, 10, p)
#' X.test[, 1] <- seq(0, 1, length.out = 10)
#' cs.pred <- predict(cs.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' cs.pred <- predict(cs.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' c.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
#'
#' # Compute a doubly robust estimate of the average treatment effect.
#' average_treatment_effect(cs.forest)
#'
#' # Compute the best linear projection on the first covariate.
#' best_linear_projection(cs.forest, X[, 1])
#'
#' # See if a causal survival forest succeeded in capturing heterogeneity by plotting
#' # the TOC and calculating a 95% CI for the AUTOC.
#' train <- sample(1:n, n / 2)
#' eval <- -train
#' train.forest <- causal_survival_forest(X[train, ], Y[train], W[train], D[train], horizon = horizon)
#' eval.forest <- causal_survival_forest(X[eval, ], Y[eval], W[eval], D[eval], horizon = horizon)
#' rate <- rank_average_treatment_effect(eval.forest,
#'                                       predict(train.forest, X[eval, ])$predictions)
#' plot(rate)
#' paste("AUTOC:", round(rate$estimate, 2), "+/", round(1.96 * rate$std.err, 2))
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

compute_psi <- function(S.hat,
                        C.hat,
                        C.Y.hat,
                        Y.hat,
                        W.centered,
                        D,
                        fY,
                        Y.index,
                        Y.grid,
                        target,
                        horizon) {
  # Compute Q(t, X) = E[f(T) | X, W, T > t]
  if (target == "RMST") {
    # Q(t, X) = E[T | X, W, T > t]
    # We can quickly compute all these t conditional expectations by updating backwards.
    # For each time point t, the conditional expectation for sample i takes the form:
    # t + Y.diff[(t + 1):grid.length] %*% S.hat[i, t:(grid.length - 1)] / S.hat[i, t]
    Y.diff <- diff(c(0, Y.grid))
    Q.hat <- matrix(NA, nrow(S.hat), ncol(S.hat))
    dot.products <- sweep(S.hat[, 1:(ncol(S.hat) - 1)], 2, Y.diff[2:ncol(S.hat)], "*")
    Q.hat[, 1] <- rowSums(dot.products)
    for (i in 2:(ncol(Q.hat) - 1)) {
      Q.hat[, i] <- Q.hat[, i - 1] - dot.products[, i - 1]
    }
    Q.hat <- Q.hat / S.hat
    Q.hat <- sweep(Q.hat, 2, Y.grid, "+") # Add back t
    Q.hat[, ncol(Q.hat)] <- max(Y.grid)
    } else {
    # Q(t, X) =  P(T > horizon | T > t)
    horizonS.index <- findInterval(horizon, Y.grid)
    Q.hat <- sweep(1 / S.hat, 1, S.hat[, horizonS.index], "*")
    Q.hat[, horizonS.index:ncol(Q.hat)] <- 1
  }

  # Pick out Q(Yi, X)
  Q.Y.hat <- Q.hat[cbind(seq_along(Y.index), Y.index)]
  numerator.one <- (D * (fY - Y.hat) + (1 - D) * (Q.Y.hat - Y.hat)) * W.centered / C.Y.hat

  # The conditional hazard function differential -d log(C.hat(t, x, w))
  # This simple forward difference approximation works reasonably well.
  # (note the "/dt" term is not needed as it cancels out in the lambda.C.hat / C.hat integral)
  log.surv.C <- -log(cbind(1, C.hat))
  dlambda.C.hat <- log.surv.C[, 2:(ncol(C.hat) + 1)] - log.surv.C[, 1:ncol(C.hat)]

  integrand <- dlambda.C.hat / C.hat * (Q.hat - Y.hat)
  numerator.two <- rep(0, length(Y.index))
  for (sample in seq_along(Y.index)) {
    Yi.index <- Y.index[sample]
    numerator.two[sample] <- sum(integrand[sample, seq_len(Yi.index)]) * W.centered[sample]
  }

  numerator <- numerator.one - numerator.two
  denominator <- W.centered^2 # denominator simplifies to this.

  list(numerator = numerator, denominator = denominator, C.Y.hat = C.Y.hat)
}
