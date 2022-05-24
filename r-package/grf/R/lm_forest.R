#' LM Forest
#'
#' Trains a linear model forest that can be used to estimate
#' \eqn{h_k(x)}, k = 1..K at X = x in the the conditional linear model
#' \eqn{Y = c(x) + h_1(x)W_1 + ... + h_K(x)W_K},
#' where Y is a (potentially vector-valued) response and
#' W a set of regressors.
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome (must be a numeric vector or matrix [one column per outcome] with no NAs).
#'  Multiple outcomes should be on the same scale.
#' @param W The conditional regressors (must be a vector or matrix with no NAs).
#' @param Y.hat Estimates of the conditional means E[Y | Xi].
#'              If Y.hat = NULL, these are estimated using
#'              a separate multi-task regression forest. Default is NULL.
#' @param W.hat Estimates of the conditional means E[Wk | Xi].
#'              If W.hat = NULL, these are estimated using
#'              a separate multi-task regression forest. Default is NULL.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param sample.weights Weights given to each sample in estimation.
#'                       If NULL, each observation receives the same weight.
#'                       Default is NULL.
#' @param gradient.weights Weights given to each coefficient h_k(x) when targeting heterogeneity
#'  in the estimates. These enter the GRF algorithm through the split criterion \eqn{\Delta}:
#'  the k-th coordinate of this is \eqn{\Delta_k} * gradient.weights[k].
#'  If NULL, each coefficient is given the same weight.
#'  Default is NULL.
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
#' @param stabilize.splits Whether or not Wk should be taken into account when
#' determining the imbalance of a split. It is an exact extension of the single-arm constraints (detailed
#' in the causal forest algorithm reference) to multiple arms, where the constraints apply to each regressor Wk.
#' Default is FALSE.
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2. (Confidence intervals are
#'                      currently only supported for univariate outcomes Y).
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained lm forest object.
#'
#' @references Athey, Susan, Julie Tibshirani, and Stefan Wager. "Generalized Random Forests".
#'  Annals of Statistics, 47(2), 2019.
#' @references Zeileis, Achim, Torsten Hothorn, and Kurt Hornik. "Model-based Recursive Partitioning."
#'  Journal of Computational and Graphical Statistics 17(2), 2008.
#'
#' @examples
#' \donttest{
#' if (require("rdd", quietly = TRUE)) {
#' # Train a LM Forest to estimate CATEs in a regression discontinuity design.
#' # Simulate a simple example with a heterogeneous jump in the CEF.
#' n <- 2000
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Z <- runif(n, -4, 4)
#' cutoff <- 0
#' W <- as.numeric(Z >= cutoff)
#' tau <- pmax(0.5 * X[, 1], 0)
#' Y <- tau * W  + 1 / (1 + exp(2 * Z)) + 0.2 * rnorm(n)
#'
#' # Compute the Imbens-Kalyanaraman MSE-optimal bandwidth for a local linear regression.
#' bandwidth <- IKbandwidth(Z, Y, cutoff)
#' # Compute kernel weights for a triangular kernel.
#' sample.weights <- kernelwts(Z, cutoff, bandwidth, "triangular")
#'
#' # Alternatively, specify bandwith and triangular kernel weights without using the `rdd` package.
#' # bandwidth <- # user can hand-specify this.
#' # dist <- abs((Z - cutoff) / bandwidth)
#' # sample.weights <- (1 - dist) * (dist <= 1) / bandwidth
#'
#' # Estimate a local linear regression with the running variable Z conditional on covariates X = x:
#' # Y = c(x) + tau(x) W + b(x) Z.
#' # Specify gradient.weights = c(1, 0) to target heterogeneity in the RDD coefficient tau(x).
#' # Also, fit forest on subset with non-zero weights for faster estimation.
#' subset <- sample.weights > 0
#' lmf <- lm_forest(X[subset, ], Y[subset], cbind(W, Z)[subset, ],
#'                  sample.weights = sample.weights[subset], gradient.weights = c(1, 0))
#' tau.hat <- predict(lmf)$predictions[, 1, ]
#'
#' # Plot estimated tau(x) vs simulated ground truth.
#' plot(X[subset, 1], tau.hat)
#' points(X[subset, 1], tau[subset], col = "red", cex = 0.1)
#' }
#' }
#'
#' @export
lm_forest <- function(X, Y, W,
                      Y.hat = NULL,
                      W.hat = NULL,
                      num.trees = 2000,
                      sample.weights = NULL,
                      gradient.weights = NULL,
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
                      stabilize.splits = FALSE,
                      ci.group.size = 2,
                      compute.oob.predictions = TRUE,
                      num.threads = NULL,
                      seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X, allow.matrix = TRUE)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)
  if (NROW(W) != nrow(X)) {
    stop("W has incorrect dimensions.")
  }
  if (!(is.numeric(W) || is.logical(W)) || anyNA(W)) {
    stop("W should be numeric with no missing values.")
  }
  if (any(apply(as.matrix(W), 2, sd) == 0)) {
    warning("W contains column(s) with zero variation (a local intercept is always included).", immediate. = TRUE)
  }
  if (is.null(gradient.weights)) {
    gradient.weights <- rep(1, NCOL(W) * NCOL(Y))
  } else {
    if (!is.numeric(gradient.weights) || length(gradient.weights) != NCOL(W)) {
      stop("gradient.weights has incorrect length (should be a numeric vector of length ncol(W)).")
    }
    gradient.weights <- rep(gradient.weights, NCOL(Y))
  }

  args.orthog <- list(X = X,
                      num.trees = max(50, num.trees / 4),
                      sample.weights = sample.weights,
                      clusters = clusters,
                      equalize.cluster.weights = equalize.cluster.weights,
                      sample.fraction = sample.fraction,
                      mtry = mtry,
                      min.node.size = 5,
                      honesty = TRUE,
                      honesty.fraction = 0.5,
                      honesty.prune.leaves = honesty.prune.leaves,
                      alpha = alpha,
                      imbalance.penalty = imbalance.penalty,
                      num.threads = num.threads,
                      seed = seed)

  if (is.null(Y.hat)) {
    forest.Y <- do.call(multi_regression_forest, c(Y = list(Y), args.orthog))
    Y.hat <- predict(forest.Y)$predictions
  } else if (is.numeric(Y.hat) && length(Y.hat) == NCOL(Y)) {
    Y.hat <- matrix(Y.hat, nrow = NROW(Y), ncol = NCOL(Y), byrow = TRUE)
  } else if (!(is.matrix(Y.hat) || is.data.frame(Y.hat) || is.vector(Y.hat))) {
      stop("Y.hat should be a matrix of E[Y | Xi] estimates.")
  } else {
    Y.hat <- as.matrix(Y.hat)
    if (NROW(Y.hat) != nrow(X) || NCOL(Y.hat) != NCOL(Y)) {
      stop("Y.hat has incorrect dimensions.")
    }
  }

  if (is.null(W.hat)) {
    forest.W <- do.call(multi_regression_forest, c(Y = list(W), args.orthog))
    W.hat <- predict(forest.W)$predictions
  } else if (is.numeric(W.hat) && length(W.hat) == NCOL(W)) {
    W.hat <- matrix(W.hat, nrow = NROW(W), ncol = NCOL(W), byrow = TRUE)
  } else if (!(is.matrix(W.hat) || is.data.frame(W.hat) || is.vector(W.hat))) {
      stop("W.hat should be a matrix of E[W | Xi] estimates.")
  } else {
    W.hat <- as.matrix(W.hat)
    if (NROW(W.hat) != nrow(X) || NCOL(W.hat) != NCOL(W)) {
      stop("W.hat has incorrect dimensions.")
    }
  }

  Y.centered <- Y - Y.hat
  W.centered <- W - W.hat
  data <- create_train_matrices(X,
                                outcome = Y.centered,
                                treatment = W.centered,
                                sample.weights = sample.weights)
  args <- list(num.trees = num.trees,
               clusters = clusters,
               samples.per.cluster = samples.per.cluster,
               gradient.weights = gradient.weights,
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

  forest <- do.call.rcpp(multi_causal_train, c(data, args))
  class(forest) <- c("lm_forest", "grf")
  forest[["seed"]] <- seed
  forest[["ci.group.size"]] <- ci.group.size
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["Y.hat"]] <- Y.hat
  forest[["W.hat"]] <- W.hat
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["sample.weights"]] <- sample.weights
  forest[["has.missing.values"]] <- has.missing.values

  forest
}

#' Predict with a lm forest
#'
#' Gets estimates of \eqn{h_k(x)}, k = 1..K in the conditionally linear model
#' \eqn{Y = c(x) + h_1(x)W_1 + ... + h_K(x)W_K}, for a target sample X = x.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{h_k}(x) are desired
#'                          (for confidence intervals). This option is currently
#'                          only supported for univariate outcomes Y.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with elements `predictions`: a 3d array of dimension [num.samples, K, M] with
#' predictions for regressor W, for each outcome 1,..,M (singleton dimensions in this array can
#' be dropped by passing the `drop` argument to `[`, or with the shorthand `$predictions[,,]`),
#'  and optionally `variance.estimates`: a matrix with K columns with variance estimates.
#'
#' @examples
#' \donttest{
#' if (require("rdd", quietly = TRUE)) {
#' # Train a LM Forest to estimate CATEs in a regression discontinuity design.
#' # Simulate a simple example with a heterogeneous jump in the CEF.
#' n <- 2000
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Z <- runif(n, -4, 4)
#' cutoff <- 0
#' W <- as.numeric(Z >= cutoff)
#' tau <- pmax(0.5 * X[, 1], 0)
#' Y <- tau * W  + 1 / (1 + exp(2 * Z)) + 0.2 * rnorm(n)
#'
#' # Compute the Imbens-Kalyanaraman MSE-optimal bandwidth for a local linear regression.
#' bandwidth <- IKbandwidth(Z, Y, cutoff)
#' # Compute kernel weights for a triangular kernel.
#' sample.weights <- kernelwts(Z, cutoff, bandwidth, "triangular")
#'
#' # Alternatively, specify bandwith and triangular kernel weights without using the `rdd` package.
#' # bandwidth <- # user can hand-specify this.
#' # dist <- abs((Z - cutoff) / bandwidth)
#' # sample.weights <- (1 - dist) * (dist <= 1) / bandwidth
#'
#' # Estimate a local linear regression with the running variable Z conditional on covariates X = x:
#' # Y = c(x) + tau(x) W + b(x) Z.
#' # Specify gradient.weights = c(1, 0) to target heterogeneity in the RDD coefficient tau(x).
#' # Also, fit forest on subset with non-zero weights for faster estimation.
#' subset <- sample.weights > 0
#' lmf <- lm_forest(X[subset, ], Y[subset], cbind(W, Z)[subset, ],
#'                  sample.weights = sample.weights[subset], gradient.weights = c(1, 0))
#' tau.hat <- predict(lmf)$predictions[, 1, ]
#'
#' # Plot estimated tau(x) vs simulated ground truth.
#' plot(X[subset, 1], tau.hat)
#' points(X[subset, 1], tau[subset], col = "red", cex = 0.1)
#' }
#' }
#'
#' @method predict lm_forest
#' @export
predict.lm_forest <- function(object,
                              newdata = NULL,
                              num.threads = NULL,
                              estimate.variance = FALSE,
                              ...) {
  if (estimate.variance && NCOL(object[["Y.orig"]]) > 1) {
    stop("Pointwise variance estimates are only supported for one outcome.")
  }
  num.outcomes <- NCOL(object[["Y.orig"]])
  outcome.names <- if (is.null(colnames(object[["Y.orig"]]))) {
    paste("Y", 1:NCOL(object[["Y.orig"]]), sep = ".")
  } else {
    make.names(colnames(object[["Y.orig"]]), unique = TRUE)
  }

  num.W <- NCOL(object[["W.orig"]])
  W.names <- if (is.null(colnames(object[["W.orig"]]))) {
    paste("h", 1:NCOL(object[["W.orig"]]), sep = ".")
  } else {
    make.names(colnames(object[["W.orig"]]), unique = TRUE)
  }

  dimnames <- list(NULL, W.names, outcome.names)
  # If possible, use pre-computed predictions.
  if (is.null(newdata) && !estimate.variance && !is.null(object$predictions)) {
    predictions <- array(object$predictions, dim = c(NROW(object$predictions), num.W, num.outcomes),
                         dimnames = dimnames)
    return(list(predictions = predictions))
  }

  num.threads <- validate_num_threads(num.threads)
  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X)

  args <- list(forest.object = forest.short,
               num.outcomes = num.outcomes,
               num.treatments = num.W,
               num.threads = num.threads,
               estimate.variance = estimate.variance)

  if (!is.null(newdata)) {
    validate_newdata(newdata, X, allow.na = TRUE)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(multi_causal_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(multi_causal_predict_oob, c(train.data, args))
  }
  predictions <- array(ret$predictions, dim = c(NROW(ret$predictions), num.W, num.outcomes),
                       dimnames = dimnames)

  list(predictions = predictions,
       variance.estimates = if (estimate.variance) ret$variance.estimates)
}
