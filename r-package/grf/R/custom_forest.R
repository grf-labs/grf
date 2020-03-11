#' Custom forest
#'
#' Trains a custom forest model.
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction. Default is 0.5.
#' @param mtry Number of variables tried for each split. Default is
#'             \eqn{\sqrt p + 20} where p is the number of variables.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#'                      Default is 5.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting). Default is TRUE.
#'  For a detailed description of honesty, honesty.fraction, honesty.prune.leaves, and recommendations for
#'  parameter tuning, see the grf
#'  \href{https://grf-labs.github.io/grf/REFERENCE.html#honesty-honesty-fraction-honesty-prune-leaves}{algorithm reference}.
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
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'  Default is NULL (ignored).
#' @param equalize.cluster.weights If FALSE, each unit is given the same weight (so that bigger
#'  clusters get more weight). If TRUE, each cluster is given equal weight in the forest. In this case,
#'  during training, each tree uses the same number of observations from each drawn cluster: If the
#'  smallest cluster has K units, then when we sample a cluster during training, we only give a random
#'  K elements of the cluster to the tree-growing procedure. When estimating average treatment effects,
#'  each observation is given weight 1/cluster size, so that the total weight of each cluster is the
#'  same.
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained regression forest object.
#'
#' @examples
#' \dontrun{
#' # Train a custom forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' c.forest <- custom_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' c.pred <- predict(c.forest, X.test)
#' }
#'
#' @export
custom_forest <- function(X, Y,
                          sample.fraction = 0.5,
                          mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                          num.trees = 2000,
                          min.node.size = 5,
                          honesty = TRUE,
                          honesty.fraction = 0.5,
                          honesty.prune.leaves = TRUE,
                          alpha = 0.05,
                          imbalance.penalty = 0.0,
                          clusters = NULL,
                          equalize.cluster.weights = FALSE,
                          compute.oob.predictions = TRUE,
                          num.threads = NULL,
                          seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X)
  Y <- validate_observations(Y, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, NULL)
  num.threads <- validate_num_threads(num.threads)

  no.split.variables <- numeric(0)

  data <- create_data_matrices(X, outcome = Y)
  ci.group.size <- 1

  forest <- custom_train(
    data$train.matrix, data$sparse.train.matrix, data$outcome.index, mtry, num.trees, min.node.size,
    sample.fraction, honesty, honesty.fraction, honesty.prune.leaves, ci.group.size, alpha,
    imbalance.penalty, clusters, samples.per.cluster, num.threads, compute.oob.predictions, seed
  )

  class(forest) <- c("custom_forest", "grf")
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["has.missing.values"]] <- has.missing.values
  forest
}

#' Predict with a custom forest.
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
#' @return Vector of predictions.
#'
#' @examples
#' \dontrun{
#' # Train a custom forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' c.forest <- custom_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' c.pred <- predict(c.forest, X.test)
#' }
#'
#' @method predict custom_forest
#' @export
predict.custom_forest <- function(object, newdata = NULL, num.threads = NULL, ...) {
  forest.short <- object[-which(names(object) == "X.orig")]

  X <- object[["X.orig"]]
  train.data <- create_data_matrices(X, object[["Y.orig"]])

  num.threads <- validate_num_threads(num.threads)

  if (!is.null(newdata)) {
    validate_newdata(newdata, X)
    data <- create_data_matrices(newdata)
    custom_predict(
      forest.short, train.data$train.matrix, train.data$sparse.train.matrix, train.data$outcome.index,
      data$train.matrix, data$sparse.train.matrix, num.threads
    )
  } else {
    custom_predict_oob(forest.short, train.data$train.matrix, train.data$sparse.train.matrix, train.data$outcome.index, num.threads)
  }
}
