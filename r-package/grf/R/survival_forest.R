#' Survival forest
#'
#' Trains a forest for right-censored surival data that can be used to estimate the
#' conditional survival function S(t, x) = P[T > t | X = x]
#'
#' @param X The covariates.
#' @param Y The event time (must be non-negative).
#' @param D The event type (0: censored, 1: failure).
#' @param failure.times A vector of event times to fit the survival curve at. If NULL, then all the observed
#'  failure times are used. This speeds up forest estimation by constraining the event grid. Observed event
#'  times are rounded down to the last sorted occurance less than or equal to the specified failure time.
#'  The time points should be in increasing order. Default is NULL.
#' @param num.trees Number of trees grown in the forest. Default is 1000.
#' @param sample.weights Weights given to an observation in prediction.
#'                       If NULL, each observation is given the same weight. Default is NULL.
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
#'                      Default is 15.
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
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. The number of failures in
#'  each child has to be at least one or `alpha` times the number of samples in the parent node. Default is 0.05.
#'  (On data with very low event rate the default value may be too high for the forest to split
#'  and lowering it may be beneficial).
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param prediction.type The type of estimate of the survival function, choices are "Kaplan-Meier" or "Nelson-Aalen".
#' Only relevant if `compute.oob.predictions` is TRUE. Default is "Kaplan-Meier".
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained survival_forest forest object.
#'
#' @references Cui, Yifan, Michael R. Kosorok, Erik Sverdrup, Stefan Wager, and Ruoqing Zhu.
#'  "Estimating Heterogeneous Treatment Effects with Right-Censored Data via Causal Survival Forests."
#'  arXiv preprint arXiv:2001.09887, 2020.
#' @references Ishwaran, Hemant, Udaya B. Kogalur, Eugene H. Blackstone, and Michael S. Lauer.
#'   "Random survival forests." The Annals of Applied Statistics 2.3 (2008): 841-860.
#'
#' @examples
#' \donttest{
#' # Train a standard survival forest.
#' n <- 2000
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' failure.time <- exp(0.5 * X[, 1]) * rexp(n)
#' censor.time <- 2 * rexp(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' s.forest <- survival_forest(X, Y, D)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 3, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 3)
#' s.pred <- predict(s.forest, X.test)
#'
#' # Plot the survival curve.
#' plot(NA, NA, xlab = "failure time", ylab = "survival function",
#'      xlim = range(s.pred$failure.times),
#'      ylim = c(0, 1))
#' for(i in 1:3) {
#'   lines(s.pred$failure.times, s.pred$predictions[i,], col = i)
#'   s.true = exp(-s.pred$failure.times / exp(0.5 * X.test[i, 1]))
#'   lines(s.pred$failure.times, s.true, col = i, lty = 2)
#' }
#'
#' # Predict on out-of-bag training samples.
#' s.pred <- predict(s.forest)
#'
#' # Plot the survival curve for the first five individuals.
#' matplot(s.pred$failure.times, t(s.pred$predictions[1:5, ]),
#'         xlab = "failure time", ylab = "survival function (OOB)",
#'         type = "l", lty = 1)
#'
#' # Train the forest on a less granular grid.
#' failure.summary <- summary(Y[D == 1])
#' events <- seq(failure.summary["Min."], failure.summary["Max."], by = 0.1)
#' s.forest.grid <- survival_forest(X, Y, D, failure.times = events)
#' s.pred.grid <- predict(s.forest.grid)
#' matpoints(s.pred.grid$failure.times, t(s.pred.grid$predictions[1:5, ]),
#'           type = "l", lty = 2)
#'
#' # Compute OOB concordance based on the mortality score in Ishwaran et al. (2008).
#' s.pred.nelson.aalen <- predict(s.forest, prediction.type = "Nelson-Aalen")
#' chf.score <- rowSums(-log(s.pred.nelson.aalen$predictions))
#' if (require("survival", quietly = TRUE)) {
#'  concordance(Surv(Y, D) ~ chf.score, reverse = TRUE)
#' }
#' }
#'
#' @export
survival_forest <- function(X, Y, D,
                            failure.times = NULL,
                            num.trees = 1000,
                            sample.weights = NULL,
                            clusters = NULL,
                            equalize.cluster.weights = FALSE,
                            sample.fraction = 0.5,
                            mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                            min.node.size = 15,
                            honesty = TRUE,
                            honesty.fraction = 0.5,
                            honesty.prune.leaves = TRUE,
                            alpha = 0.05,
                            prediction.type = c("Kaplan-Meier", "Nelson-Aalen"),
                            compute.oob.predictions = TRUE,
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  if (any(Y < 0)) {
    stop("The event times must be non-negative.")
  }
  D <- validate_observations(D, X)
  if (!all(D %in% c(0, 1))) {
    stop("The censor values can only be 0 or 1.")
  }
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)
  prediction.type <- match.arg(prediction.type)
  if (prediction.type == "Kaplan-Meier") {
    prediction.type <- 0
  } else if (prediction.type == "Nelson-Aalen") {
    prediction.type <- 1
  }

  # Relabel the times to consecutive integers such that:
  # if the event time is less than the smallest failure time: set it to 0
  # if the event time is above the latter, but less than the second smallest failure time: set it to 1
  # etc. Will range from 0 to num.failures.
  if (is.null(failure.times)) {
    failure.times <- sort(unique(Y[D == 1]))
  } else if (is.unsorted(failure.times, strictly = TRUE)) {
    stop("Argument `failure.times` should be a vector with elements in increasing order.")
  }
  Y.relabeled <- findInterval(Y, failure.times)

  data <- create_train_matrices(X, outcome = Y.relabeled, sample.weights = sample.weights, censor = D)
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
               num.failures = length(failure.times),
               prediction.type = prediction.type,
               compute.oob.predictions = compute.oob.predictions,
               num.threads = num.threads,
               seed = seed)

  forest <- do.call.rcpp(survival_train, c(data, args))
  class(forest) <- c("survival_forest", "grf")
  forest[["seed"]] <- seed
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["Y.relabeled"]] <- Y.relabeled
  forest[["D.orig"]] <- D
  forest[["sample.weights"]] <- sample.weights
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["has.missing.values"]] <- has.missing.values
  forest[["failure.times"]] <- failure.times
  forest[["prediction.type"]] <- prediction.type

  forest
}

#' Predict with a survival forest
#'
#' Gets estimates of the conditional survival function S(t, x) using a trained survival forest. The curve can be
#' estimated by Kaplan-Meier, or Nelson-Aalen.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param failure.times A vector of failure times to make predictions at. If NULL, then the
#'  failure times used for training the forest is used. The time points should be in increasing order. Default is NULL.
#' @param prediction.type The type of estimate of the survival function, choices are "Kaplan-Meier" or "Nelson-Aalen".
#'  The default is the prediction.type used to train the forest.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with elements `failure.times`: a vector of event times t for the survival curve,
#'  and `predictions`: a matrix of survival curves. Each row is the survival curve for
#'  sample X_i: predictions[i, j] = S(failure.times[j], X_i).
#'
#' @examples
#' \donttest{
#' # Train a standard survival forest.
#' n <- 2000
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' failure.time <- exp(0.5 * X[, 1]) * rexp(n)
#' censor.time <- 2 * rexp(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' s.forest <- survival_forest(X, Y, D)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 3, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 3)
#' s.pred <- predict(s.forest, X.test)
#'
#' # Plot the survival curve.
#' plot(NA, NA, xlab = "failure time", ylab = "survival function",
#'      xlim = range(s.pred$failure.times),
#'      ylim = c(0, 1))
#' for(i in 1:3) {
#'   lines(s.pred$failure.times, s.pred$predictions[i,], col = i)
#'   s.true = exp(-s.pred$failure.times / exp(0.5 * X.test[i, 1]))
#'   lines(s.pred$failure.times, s.true, col = i, lty = 2)
#' }
#'
#' # Predict on out-of-bag training samples.
#' s.pred <- predict(s.forest)
#'
#' # Plot the survival curve for the first five individuals.
#' matplot(s.pred$failure.times, t(s.pred$predictions[1:5, ]),
#'         xlab = "failure time", ylab = "survival function (OOB)",
#'         type = "l", lty = 1)
#'
#' # Train the forest on a less granular grid.
#' failure.summary <- summary(Y[D == 1])
#' events <- seq(failure.summary["Min."], failure.summary["Max."], by = 0.1)
#' s.forest.grid <- survival_forest(X, Y, D, failure.times = events)
#' s.pred.grid <- predict(s.forest.grid)
#' matpoints(s.pred.grid$failure.times, t(s.pred.grid$predictions[1:5, ]),
#'           type = "l", lty = 2)
#'
#' # Compute OOB concordance based on the mortality score in Ishwaran et al. (2008).
#' s.pred.nelson.aalen <- predict(s.forest, prediction.type = "Nelson-Aalen")
#' chf.score <- rowSums(-log(s.pred.nelson.aalen$predictions))
#' if (require("survival", quietly = TRUE)) {
#'  concordance(Surv(Y, D) ~ chf.score, reverse = TRUE)
#' }
#' }
#'
#' @method predict survival_forest
#' @export
predict.survival_forest <- function(object,
                                    newdata = NULL,
                                    failure.times = NULL,
                                    prediction.type = c("Kaplan-Meier", "Nelson-Aalen"),
                                    num.threads = NULL, ...) {
  num.threads <- validate_num_threads(num.threads)
  if (is.null(failure.times)) {
    failure.times <- object[["failure.times"]]
    Y.relabeled <- object[["Y.relabeled"]]
  } else {
    if (is.unsorted(failure.times, strictly = TRUE)) {
      stop("Argument `failure.times` should be a vector with elements in increasing order.")
    }
    Y.relabeled <- findInterval(object[["Y.orig"]], failure.times)
  }

  default.prediction.type <- length(prediction.type) == 2
  prediction.type <- match.arg(prediction.type)
  if (default.prediction.type) {
    prediction.type <- object[["prediction.type"]]
  } else if (prediction.type == "Kaplan-Meier") {
    prediction.type <- 0
  } else if (prediction.type == "Nelson-Aalen") {
    prediction.type <- 1
  }

  # If possible, use pre-computed predictions.
  failure.times.orig <- object[["failure.times"]]
  prediction.type.orig <- object[["prediction.type"]]
  if (is.null(newdata) && identical(failure.times, failure.times.orig)
      && identical(prediction.type, prediction.type.orig) && !is.null(object$predictions)) {
    return(list(predictions = object$predictions, failure.times = failure.times))
  }

  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X,
                                      outcome = Y.relabeled,
                                      censor = object[["D.orig"]],
                                      sample.weights = object[["sample.weights"]])

  args <- list(forest.object = forest.short,
               num.threads = num.threads,
               num.failures = length(failure.times),
               prediction.type = prediction.type)

  if (!is.null(newdata)) {
    validate_newdata(newdata, X, allow.na = TRUE)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(survival_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(survival_predict_oob, c(train.data, args))
  }

  list(predictions = ret[["predictions"]], failure.times = failure.times)
}
