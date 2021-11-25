#' Multi-task regression forest
#'
#' Trains a regression forest that can be used to estimate
#' the conditional mean functions mu_i(x) = E[Y_i | X = x]
#'
#' @param X The covariates used in the regression.
#' @param Y The outcomes (must be a numeric vector/matrix with no NAs).
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
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained multi regression forest object.
#'
#' @examples
#' \donttest{
#' # Train a standard regression forest.
#' n <- 500
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Y <-  X[, 1, drop = FALSE] %*% cbind(1, 2) + rnorm(n)
#' mr.forest <- multi_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' mr.pred <- predict(mr.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' mr.pred <- predict(mr.forest)
#' }
#'
#' @export
multi_regression_forest <- function(X, Y,
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
                                    compute.oob.predictions = TRUE,
                                    num.threads = NULL,
                                    seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X, allow.matrix = TRUE)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)

  data <- create_train_matrices(X, outcome = Y, sample.weights = sample.weights)
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
               compute.oob.predictions = compute.oob.predictions,
               num.threads = num.threads,
               seed = seed)

  forest <- do.call.rcpp(multi_regression_train, c(data, args))
  class(forest) <- c("multi_regression_forest", "grf")
  forest[["seed"]] <- seed
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["sample.weights"]] <- sample.weights
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["has.missing.values"]] <- has.missing.values

  forest
}

#' Predict with a multi regression forest
#'
#' Gets estimates of E[Y_i | X = x] using a trained multi regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list containing `predictions`: a matrix of predictions for each outcome.
#'
#' @examples
#' \donttest{
#' # Train a standard regression forest.
#' n <- 500
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Y <-  X[, 1, drop = FALSE] %*% cbind(1, 2) + rnorm(n)
#' mr.forest <- multi_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' mr.pred <- predict(mr.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' mr.pred <- predict(mr.forest)
#' }
#'
#' @method predict multi_regression_forest
#' @export
predict.multi_regression_forest <- function(object,
                                            newdata = NULL,
                                            num.threads = NULL,
                                            ...) {
  outcome.names <- if (is.null(colnames(object[["Y.orig"]]))) {
    paste0("Y", 1:NCOL(object[["Y.orig"]]))
  } else {
    colnames(object[["Y.orig"]])
  }
  # If possible, use pre-computed predictions.
  if (is.null(newdata) && !is.null(object$predictions)) {
    colnames(object$predictions) <- outcome.names
    return(list(predictions = object$predictions))
  }

  num.threads <- validate_num_threads(num.threads)
  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X)

  args <- list(forest.object = forest.short,
               num.threads = num.threads,
               num.outcomes = NCOL(object[["Y.orig"]]))

  if (!is.null(newdata)) {
    validate_newdata(newdata, X, allow.na = TRUE)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(multi_regression_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(multi_regression_predict_oob, c(train.data, args))
  }
  colnames(ret$predictions) <- outcome.names

  list(predictions = ret$predictions)
}
