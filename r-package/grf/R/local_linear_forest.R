#' Local Linear forest
#'
#' Trains a local linear forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
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
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 1.
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if FALSE
#'                        NULL parameters are set to defaults. Default is FALSE.
#' @param tune.num.trees The number of trees in each 'mini forest' used to fit the tuning model. Default is 10.
#' @param tune.num.reps The number of forests used to fit the tuning model. Default is 100.
#' @param tune.num.draws The number of random parameter values considered when using the model
#'                          to select the optimal parameters. Default is 1000.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained local linear forest object.
#'
#' @examples
#' \dontrun{
#' # Train a standard regression forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' forest <- ll_regression_forest(X, Y)
#' }
#'
#' @export
ll_regression_forest <- function(X, Y,
                                 num.trees = 2000,
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
                                 ci.group.size = 1,
                                 tune.parameters = "none",
                                 tune.num.trees = 10,
                                 tune.num.reps = 100,
                                 tune.num.draws = 1000,
                                 num.threads = NULL,
                                 seed = runif(1, 0, .Machine$integer.max)) {

  forest <- regression_forest(X, Y,
                              num.trees = num.trees,
                              sample.weights = NULL,
                              clusters = clusters,
                              equalize.cluster.weights = equalize.cluster.weights,
                              sample.fraction = sample.fraction,
                              mtry = mtry,
                              min.node.size = min.node.size,
                              honesty = honesty,
                              honesty.fraction = honesty.fraction,
                              honesty.prune.leaves = honesty.prune.leaves,
                              alpha = alpha,
                              imbalance.penalty = imbalance.penalty,
                              ci.group.size = ci.group.size,
                              tune.parameters = tune.parameters,
                              tune.num.trees = tune.num.trees,
                              tune.num.reps = tune.num.reps,
                              tune.num.draws = tune.num.draws,
                              compute.oob.predictions = FALSE,
                              num.threads = num.threads,
                              seed = seed)
  forest[["tuning.output"]] <- NULL
  class(forest) <- c("ll_regression_forest", "grf")

  forest
}

#' Predict with a local linear forest
#'
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param linear.correction.variables Optional subset of indexes for variables to be used in local
#'                   linear prediction. If left NULL, all variables are used.
#'                   We run a locally weighted linear regression on the included variables.
#'                   Please note that this is a beta feature still in development, and may slow down
#'                   prediction considerably. Defaults to NULL.
#' @param ll.lambda Ridge penalty for local linear predictions
#' @param ll.weight.penalty Option to standardize ridge penalty by covariance (TRUE),
#'                            or penalize all covariates equally (FALSE). Defaults to FALSE.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A vector of predictions.
#'
#' @examples
#' \dontrun{
#' # Train the forest.
#' n <- 50
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' forest <- ll_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' predictions <- predict(forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' predictions.oob <- predict(forest)
#' }
#'
#' @method predict ll_regression_forest
#' @export
predict.ll_regression_forest <- function(object, newdata = NULL,
                                         linear.correction.variables = NULL,
                                         ll.lambda = NULL,
                                         ll.weight.penalty = FALSE,
                                         num.threads = NULL,
                                         estimate.variance = FALSE,
                                         ...) {
  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  if (is.null(linear.correction.variables)) {
    linear.correction.variables <- 1:ncol(X)
  }
  # Validate and account for C++ indexing
  linear.correction.variables <- validate_ll_vars(linear.correction.variables, ncol(X))

  if (is.null(ll.lambda)) {
    ll.regularization.path <- tune_ll_regression_forest(
      object, linear.correction.variables,
      ll.weight.penalty, num.threads
    )
    ll.lambda <- ll.regularization.path$lambda.min
  } else {
    ll.lambda <- validate_ll_lambda(ll.lambda)
  }

  num.threads <- validate_num_threads(num.threads)

  # Subtract 1 to account for C++ indexing
  linear.correction.variables <- linear.correction.variables - 1

  train.data <- create_data_matrices(X, outcome = object[["Y.orig"]])

  if (!is.null(newdata)) {
    validate_newdata(newdata, X)
    data <- create_data_matrices(newdata)
    ret <- ll_regression_predict(
      forest.short, train.data$train.matrix, train.data$sparse.train.matrix, train.data$outcome.index,
      data$train.matrix, data$sparse.train.matrix,
      ll.lambda, ll.weight.penalty, linear.correction.variables, num.threads, estimate.variance
    )
  } else {
    ret <- ll_regression_predict_oob(
      forest.short, train.data$train.matrix, train.data$sparse.train.matrix, train.data$outcome.index,
      ll.lambda, ll.weight.penalty, linear.correction.variables, num.threads, estimate.variance
    )
  }

  ret[["ll.lambda"]] <- ll.lambda

  # Convert list to data frame.
  empty <- sapply(ret, function(elem) length(elem) == 0)
  do.call(cbind.data.frame, ret[!empty])
}
