#' Quantile forest
#'
#' Trains a regression forest that can be used to estimate
#' quantiles of the conditional distribution of Y given X = x.
#'
#' @param X The covariates used in the quantile regression.
#' @param Y The outcome.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param quantiles Vector of quantiles used to calibrate the forest. Default is (0.1, 0.5, 0.9).
#' @param regression.splitting Whether to use regression splits when growing trees instead
#'                             of specialized splits based on the quantiles (the default).
#'                             Setting this flag to true corresponds to the approach to
#'                             quantile forests from Meinshausen (2006). Default is FALSE.
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
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained quantile forest object.
#'
#' @examples
#' \dontrun{
#' # Generate data.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' Y <- X[, 1] * rnorm(n)
#'
#' # Train a quantile forest.
#' q.forest <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9))
#'
#' # Make predictions.
#' q.hat <- predict(q.forest, X.test)
#'
#' # Make predictions for different quantiles than those used in training.
#' q.hat <- predict(q.forest, X.test, quantiles = c(0.1, 0.9))
#'
#' # Train a quantile forest using regression splitting instead of quantile-based
#' # splits, emulating the approach in Meinshausen (2006).
#' meins.forest <- quantile_forest(X, Y, regression.splitting = TRUE)
#'
#' # Make predictions for the desired quantiles.
#' q.hat <- predict(meins.forest, X.test, quantiles = c(0.1, 0.5, 0.9))
#' }
#'
#' @export
quantile_forest <- function(X, Y,
                            num.trees = 2000,
                            quantiles = c(0.1, 0.5, 0.9),
                            regression.splitting = FALSE,
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
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max)) {
  if (!is.numeric(quantiles) | length(quantiles) < 1) {
    stop("Error: Must provide numeric quantiles")
  } else if (min(quantiles) <= 0 | max(quantiles) >= 1) {
    stop("Error: Quantiles must be in (0, 1)")
  }

  has.missing.values <- validate_X(X, allow.na = TRUE)
  Y <- validate_observations(Y, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, NULL)
  num.threads <- validate_num_threads(num.threads)

  data <- create_data_matrices(X, outcome = Y)
  ci.group.size <- 1

  forest <- quantile_train(
    quantiles, regression.splitting, data$train.matrix, data$sparse.train.matrix, data$outcome.index, mtry,
    num.trees, min.node.size, sample.fraction, honesty, honesty.fraction, honesty.prune.leaves,
    ci.group.size, alpha, imbalance.penalty, clusters, samples.per.cluster, num.threads, seed
  )

  class(forest) <- c("quantile_forest", "grf")
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["has.missing.values"]] <- has.missing.values
  forest
}

#' Predict with a quantile forest
#'
#' Gets estimates of the conditional quantiles of Y given X using a trained forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param quantiles Vector of quantiles at which estimates are required.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Predictions at each test point for each desired quantile.
#'
#' @examples
#' \dontrun{
#' # Train a quantile forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' q.forest <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9))
#'
#' # Predict on out-of-bag training samples.
#' q.pred <- predict(q.forest)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' q.pred <- predict(q.forest, X.test)
#' }
#'
#' @method predict quantile_forest
#' @export
predict.quantile_forest <- function(object,
                                    newdata = NULL,
                                    quantiles = c(0.1, 0.5, 0.9),
                                    num.threads = NULL, ...) {
  if (!is.numeric(quantiles) | length(quantiles) < 1) {
    stop("Error: Must provide numeric quantiles")
  } else if (min(quantiles) <= 0 | max(quantiles) >= 1) {
    stop("Error: Quantiles must be in (0, 1)")
  }

  num.threads <- validate_num_threads(num.threads)

  forest.short <- object[-which(names(object) == "X.orig")]

  X <- object[["X.orig"]]
  train.data <- create_data_matrices(X, outcome = object[["Y.orig"]])

  if (!is.null(newdata)) {
    validate_newdata(newdata, object$X.orig, allow.na = TRUE)
    data <- create_data_matrices(newdata)
    quantile_predict(
      forest.short, quantiles, train.data$train.matrix, train.data$sparse.train.matrix,
      train.data$outcome.index, data$train.matrix, data$sparse.train.matrix, num.threads
    )
  } else {
    quantile_predict_oob(
      forest.short, quantiles, train.data$train.matrix, train.data$sparse.train.matrix,
      train.data$outcome.index, num.threads
    )
  }
}
