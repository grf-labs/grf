#' Local Linear forest
#'
#' Trains a local linear forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]
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
#' @param honesty Whether or not honest splitting (i.e., sub-sample splitting) should be used. Default is TRUE.
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. When using the defaults (honesty = TRUE and
#'                         honesty.fraction = NULL), half of the data will be used for determining splits.
#'                         Default is 0.5.
#' @param prune.empty.leaves (experimental) If true, prunes the estimation sample tree such that no leaves
#'  are empty. If false, keep the same tree as determined in the splits sample (if an empty leave is encountered, that
#'  tree is skipped and does not contribute to the estimate). Setting this to false may improve performance on
#'  small/marginally powered data, but requires more trees. Only applies if honesty is enabled. Default is TRUE.
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 1.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. Default is 0.05.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized. Default is 0.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'                 Default is NULL (ignored).
#' @param samples.per.cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster when training a tree. If NULL, we set samples.per.cluster to the size
#'                            of the smallest cluster. If some clusters are smaller than samples.per.cluster,
#'                            the whole cluster is used every time the cluster is drawn. Note that
#'                            clusters with less than samples.per.cluster observations get relatively
#'                            smaller weight than others in training the forest, i.e., the contribution
#'                            of a given cluster to the final forest scales with the minimum of
#'                            the number of observations in the cluster and samples.per.cluster. Default is NULL.
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if false
#'                        NULL parameters are set to defaults. Default is FALSE.
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model. Default is 10.
#' @param num.fit.reps The number of forests used to fit the tuning model. Default is 100.
#' @param num.optimize.reps The number of random parameter values considered when using the model
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
                                 sample.fraction = 0.5,
                                 mtry = NULL,
                                 num.trees = 2000,
                                 min.node.size = NULL,
                                 honesty = TRUE,
                                 honesty.fraction = NULL,
                                 prune.empty.leaves = NULL,
                                 ci.group.size = 1,
                                 alpha = NULL,
                                 imbalance.penalty = NULL,
                                 clusters = NULL,
                                 samples.per.cluster = NULL,
                                 tune.parameters = FALSE,
                                 num.fit.trees = 10,
                                 num.fit.reps = 100,
                                 num.optimize.reps = 1000,
                                 num.threads = NULL,
                                 seed = NULL) {
  validate_X(X)
  Y <- validate_observations(Y, X)

  num.threads <- validate_num_threads(num.threads)
  seed <- validate_seed(seed)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_samples_per_cluster(samples.per.cluster, clusters)

  if (tune.parameters) {
    tuning.output <- tune_regression_forest(X, Y,
      num.fit.trees = num.fit.trees,
      num.fit.reps = num.fit.reps,
      num.optimize.reps = num.optimize.reps,
      min.node.size = min.node.size,
      sample.fraction = sample.fraction,
      mtry = mtry,
      alpha = alpha,
      imbalance.penalty = imbalance.penalty,
      num.threads = num.threads,
      honesty = honesty,
      honesty.fraction = honesty.fraction,
      prune.empty.leaves = prune.empty.leaves,
      seed = seed,
      clusters = clusters,
      samples.per.cluster = samples.per.cluster
    )
    tunable.params <- tuning.output$params
  } else {
    tunable.params <- c(
      min.node.size = validate_min_node_size(min.node.size),
      sample.fraction = validate_sample_fraction(sample.fraction),
      mtry = validate_mtry(mtry, X),
      alpha = validate_alpha(alpha),
      imbalance.penalty = validate_imbalance_penalty(imbalance.penalty),
      honesty.fraction = validate_honesty_fraction(honesty.fraction, honesty),
      prune.empty.leaves = validate_prune_empty_leaves(prune.empty.leaves)
    )
  }

  data <- create_data_matrices(X, Y)
  outcome.index <- ncol(X) + 1
  sample.weight.index <- ncol(X) + 2
  compute.oob.predictions <- FALSE

  forest <- regression_train(
    data$default, data$sparse, outcome.index, sample.weight.index,
    FALSE,
    as.numeric(tunable.params["mtry"]),
    num.trees,
    as.numeric(tunable.params["min.node.size"]),
    as.numeric(tunable.params["sample.fraction"]),
    honesty,
    as.numeric(tunable.params["honesty.fraction"]),
    as.numeric(tunable.params["prune.empty.leaves"]),
    ci.group.size,
    as.numeric(tunable.params["alpha"]),
    as.numeric(tunable.params["imbalance.penalty"]),
    clusters,
    samples.per.cluster,
    compute.oob.predictions,
    num.threads,
    seed
  )

  class(forest) <- c("ll_regression_forest", "grf")
  forest[["ci.group.size"]] <- ci.group.size
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["clusters"]] <- clusters
  forest[["tunable.params"]] <- tunable.params
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

  train.data <- create_data_matrices(X, object[["Y.orig"]])
  outcome.index <- ncol(X) + 1

  if (!is.null(newdata)) {
    validate_newdata(newdata, X)
    data <- create_data_matrices(newdata)
    ret <- ll_regression_predict(
      forest.short, train.data$default, train.data$sparse, outcome.index,
      data$default, data$sparse,
      ll.lambda, ll.weight.penalty, linear.correction.variables, num.threads, estimate.variance
    )
  } else {
    ret <- ll_regression_predict_oob(
      forest.short, train.data$default, train.data$sparse, outcome.index,
      ll.lambda, ll.weight.penalty, linear.correction.variables, num.threads, estimate.variance
    )
  }

  ret[["ll.lambda"]] <- ll.lambda

  # Convert list to data frame.
  empty <- sapply(ret, function(elem) length(elem) == 0)
  do.call(cbind.data.frame, ret[!empty])
}
