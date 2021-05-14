#' Probability forest
#'
#' Trains a probability forest that can be used to estimate
#' the conditional class probabilities P[Y = k | X = x]
#'
#' @param X The covariates.
#' @param Y The class label (must be a factor vector with no NAs).
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param sample.weights Weights given to an observation in estimation.
#'                       If NULL, each observation is given the same weight. Default is NULL.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'  Default is NULL (ignored).
#' @param equalize.cluster.weights If FALSE, each unit is given the same weight (so that bigger
#'  clusters get more weight). If TRUE, each cluster is given equal weight in the forest. In this case,
#'  during training, each tree uses the same number of observations from each drawn cluster: If the
#'  smallest cluster has K units, then when we sample a cluster during training, we only give a random
#'  K elements of the cluster to the tree-growing procedure. When estimating average treatment effects,
#'  each observation is given weight 1/cluster size, so that the total weight of each cluster is the
#'  same.
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
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2.
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained probability forest object.
#'
#' @examples
#' \donttest{
#' # Train a probability forest.
#' p <- 5
#' n <- 2000
#' X <- matrix(rnorm(n*p), n, p)
#' prob <- 1 / (1 + exp(-X[, 1] - X[, 2]))
#' Y <- as.factor(rbinom(n, 1, prob))
#' p.forest <- probability_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 10, p)
#' X.test[, 1] <- seq(-1.5, 1.5, length.out = 10)
#' p.hat <- predict(p.forest, X.test, estimate.variance = TRUE)
#'
#' # Plot the estimated success probabilities with 95 % confidence bands.
#' prob.test <- 1 / (1 + exp(-X.test[, 1] - X.test[, 2]))
#' p.true <- cbind(`0` = 1 - prob.test, `1` = prob.test)
#' plot(X.test[, 1], p.true[, "1"], col = 'red', ylim = c(0, 1))
#' points(X.test[, 1], p.hat$predictions[, "1"], pch = 16)
#' lines(X.test[, 1], (p.hat$predictions + 2 * sqrt(p.hat$variance.estimates))[, "1"])
#' lines(X.test[, 1], (p.hat$predictions - 2 * sqrt(p.hat$variance.estimates))[, "1"])
#'
#' # Predict on out-of-bag training samples.
#' p.hat <- predict(p.forest)
#' }
#'
#' @export
probability_forest <- function(X, Y,
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
                               imbalance.penalty = 0.0,
                               ci.group.size = 2,
                               compute.oob.predictions = TRUE,
                               num.threads = NULL,
                               seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)
  if (length(Y) != nrow(X)) {
    stop("length of observations Y does not equal nrow(X).")
  }
  if (anyNA(Y)) {
    stop("The vector of observations contains at least one NA.")
  }
  if (!is.factor(Y)) {
    stop("The class labels must be a factor vector.")
  }
  Y.relabeled <- as.numeric(Y) - 1 # convert to integers between 0 and num_classes
  class.names <- levels(Y)
  num.classes <- length(class.names)

  data <- create_train_matrices(X, outcome = Y.relabeled, sample.weights = sample.weights)
  args <- list(num.classes = num.classes,
               num.trees = num.trees,
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
               ci.group.size = ci.group.size,
               compute.oob.predictions = compute.oob.predictions,
               num.threads = num.threads,
               seed = seed)

  forest <- do.call.rcpp(probability_train, c(data, args))
  class(forest) <- c("probability_forest", "grf")
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["Y.relabeled"]] <- Y.relabeled
  forest[["sample.weights"]] <- sample.weights
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["has.missing.values"]] <- has.missing.values
  forest[["num.classes"]] <- num.classes
  forest[["class.names"]] <- class.names

  forest
}

#' Predict with a probability forest
#'
#' Gets estimates of P[Y = k | X = x] using a trained forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for P[Y = k | X] are desired (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with attributes `predictions`: a matrix of predictions for each class, and optionally
#'  the attribute `variance.estimates`: a matrix of variance estimates for each class.
#'
#' @examples
#' \donttest{
#' # Train a probability forest.
#' p <- 5
#' n <- 2000
#' X <- matrix(rnorm(n*p), n, p)
#' prob <- 1 / (1 + exp(-X[, 1] - X[, 2]))
#' Y <- as.factor(rbinom(n, 1, prob))
#' p.forest <- probability_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 10, p)
#' X.test[, 1] <- seq(-1.5, 1.5, length.out = 10)
#' p.hat <- predict(p.forest, X.test, estimate.variance = TRUE)
#'
#' # Plot the estimated success probabilities with 95 % confidence bands.
#' prob.test <- 1 / (1 + exp(-X.test[, 1] - X.test[, 2]))
#' p.true <- cbind(`0` = 1 - prob.test, `1` = prob.test)
#' plot(X.test[, 1], p.true[, "1"], col = 'red', ylim = c(0, 1))
#' points(X.test[, 1], p.hat$predictions[, "1"], pch = 16)
#' lines(X.test[, 1], (p.hat$predictions + 2 * sqrt(p.hat$variance.estimates))[, "1"])
#' lines(X.test[, 1], (p.hat$predictions - 2 * sqrt(p.hat$variance.estimates))[, "1"])
#'
#' # Predict on out-of-bag training samples.
#' p.hat <- predict(p.forest)
#' }
#'
#' @method predict probability_forest
#' @export
predict.probability_forest <- function(object,
                                       newdata = NULL,
                                       num.threads = NULL,
                                       estimate.variance = FALSE, ...) {
  class.names <- object$class.names
  # If possible, use pre-computed predictions.
  if (is.null(newdata) && !estimate.variance && !is.null(object$predictions)) {
    colnames(object$predictions) <- class.names
    return(list(predictions = object$predictions))
  }

  num.threads <- validate_num_threads(num.threads)
  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X, outcome = object[["Y.relabeled"]])

  args <- list(forest.object = forest.short,
               num.classes = object[["num.classes"]],
               num.threads = num.threads,
               estimate.variance = estimate.variance)

  if (!is.null(newdata)) {
    validate_newdata(newdata, object$X.orig, allow.na = TRUE)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(probability_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(probability_predict_oob, c(train.data, args))
  }
  colnames(ret$predictions) <- class.names

  list(predictions = ret$predictions,
       variance.estimates = if (estimate.variance) ret$variance.estimates)
}
