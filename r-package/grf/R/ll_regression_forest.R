#' Local linear forest
#'
#' Trains a local linear forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param enable.ll.split (experimental) Optional choice to make forest splits based on ridge residuals as opposed to
#'                        standard CART splits. Defaults to FALSE.
#' @param ll.split.weight.penalty If using local linear splits, user can specify whether or not to use a
#'                                covariance ridge penalty, analogously to the prediction case. Defaults to FALSE.
#' @param ll.split.lambda Ridge penalty for splitting. Defaults to 0.1.
#' @param ll.split.variables Linear correction variables for splitting. Defaults to all variables.
#' @param ll.split.cutoff Enables the option to use regression coefficients from the full dataset for LL splitting
#'                        once leaves get sufficiently small. Leaf size after which we use the overall beta.
#'                        Defaults to the square root of the number of samples. If desired, users can enforce no
#'                        regulation (i.e., using the leaf betas at each step) by setting this parameter to zero.
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
#'  same. Default is FALSE.
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
#'                      be at least 2. Default is 1.
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if FALSE
#'                        NULL parameters are set to defaults. Default is FALSE. Currently, local linear tuning
#'                        is based on regression forest fit, and is only supported for `enable.ll.split = FALSE`.
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
#' \donttest{
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
                                enable.ll.split = FALSE,
                                ll.split.weight.penalty = FALSE,
                                ll.split.lambda = 0.1,
                                ll.split.variables = NULL,
                                ll.split.cutoff = NULL,
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
                                ci.group.size = 2,
                                tune.parameters = "none",
                                tune.num.trees = 50,
                                tune.num.reps = 100,
                                tune.num.draws = 1000,
                                num.threads = NULL,
                                seed = runif(1, 0, .Machine$integer.max)) {

  has.missing.values <- validate_X(X)
  Y <- validate_observations(Y, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, NULL)
  num.threads <- validate_num_threads(num.threads)

  ll.split.variables <- validate_ll_vars(ll.split.variables, ncol(X))
  ll.split.lambda <- validate_ll_lambda(ll.split.lambda)
  ll.split.cutoff <- validate_ll_cutoff(ll.split.cutoff, nrow(X))

  all.tunable.params <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
                          "honesty.prune.leaves", "alpha", "imbalance.penalty")
  default.parameters <- list(sample.fraction = 0.5,
                             mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                             min.node.size = 5,
                             honesty.fraction = 0.5,
                             honesty.prune.leaves = TRUE,
                             alpha = 0.05,
                             imbalance.penalty = 0)

  # The ll_regression train wrapper signature does not contain sample weights, which is why sample.weights=FALSE,
  # whereas the regression_train wrapper signature does, and specifying sample.weights=NULL disables them.
  data <- create_train_matrices(X, outcome = Y, sample.weights = if (enable.ll.split) FALSE else NULL)

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
               ci.group.size = ci.group.size,
               num.threads = num.threads,
               seed = seed)
  if (enable.ll.split && ll.split.cutoff > 0) {
    # find overall beta
    J <- diag(ncol(X) + 1)
    J[1,1] <- 0
    D <- cbind(1, X)
    overall.beta <- solve(t(D) %*% D + ll.split.lambda * J) %*% t(D) %*% Y

    # update arguments with LLF parameters
    args <- c(args, list(ll.split.weight.penalty = ll.split.weight.penalty,
                         ll.split.lambda = ll.split.lambda,
                         ll.split.variables = ll.split.variables,
                         ll.split.cutoff = ll.split.cutoff,
                         overall.beta = overall.beta))
  } else if (enable.ll.split) {
    # update arguments with LLF parameters
    args <- c(args, list(ll.split.weight.penalty = ll.split.weight.penalty,
                         ll.split.lambda = ll.split.lambda,
                         ll.split.variables = ll.split.variables,
                         ll.split.cutoff = ll.split.cutoff,
                         overall.beta = vector(mode = "numeric", length = 0)))
  } else {
    args <- c(args, compute.oob.predictions = FALSE)
  }

  tuning.output <- NULL
  if (!identical(tune.parameters, "none")) {
    if (enable.ll.split) {
      stop("Tuning is currently only supported when enable.ll.split = FALSE.")
    }
    if (identical(tune.parameters, "all")) {
      tune.parameters <- all.tunable.params
    } else {
      tune.parameters <- unique(match.arg(tune.parameters, all.tunable.params, several.ok = TRUE))
    }
    if (!honesty) {
      tune.parameters <- tune.parameters[!grepl("honesty", tune.parameters)]
    }
    tune.parameters.defaults <- default.parameters[tune.parameters]
    tuning.output <- tune_forest(data = data,
                                 nrow.X = nrow(X),
                                 ncol.X = ncol(X),
                                 args = args,
                                 tune.parameters = tune.parameters,
                                 tune.parameters.defaults = tune.parameters.defaults,
                                 tune.num.trees = tune.num.trees,
                                 tune.num.reps = tune.num.reps,
                                 tune.num.draws = tune.num.draws,
                                 train = regression_train)

    args <- utils::modifyList(args, as.list(tuning.output[["params"]]))
  }

  if (enable.ll.split) {
    forest <- do.call.rcpp(ll_regression_train, c(data, args))
  } else {
    forest <- do.call.rcpp(regression_train, c(data, args))
  }

  class(forest) <- c("ll_regression_forest", "grf")
  forest[["seed"]] <- seed
  forest[["ci.group.size"]] <- ci.group.size
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["tunable.params"]] <- args[all.tunable.params]
  forest[["tuning.output"]] <- tuning.output
  forest[["has.missing.values"]] <- has.missing.values

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
#' @param ll.lambda Ridge penalty for local linear predictions. Defaults to NULL and will be cross-validated.
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
#' \donttest{
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
  num.threads <- validate_num_threads(num.threads)
  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X, outcome = object[["Y.orig"]])

  linear.correction.variables <- validate_ll_vars(linear.correction.variables, ncol(X))
  if (is.null(ll.lambda)) {
    ll.regularization.path <- tune_ll_regression_forest(
      object, linear.correction.variables,
      ll.weight.penalty, num.threads)
    ll.lambda <- ll.regularization.path$lambda.min
  } else {
    ll.lambda <- validate_ll_lambda(ll.lambda)
  }
  # Subtract 1 to account for C++ indexing
  linear.correction.variables <- linear.correction.variables - 1
  args <- list(forest.object = forest.short,
               num.threads = num.threads,
               estimate.variance = estimate.variance,
               ll.lambda = ll.lambda,
               ll.weight.penalty = ll.weight.penalty,
               linear.correction.variables = linear.correction.variables)

  if (!is.null(newdata)) {
    validate_newdata(newdata, X)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(ll_regression_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(ll_regression_predict_oob, c(train.data, args))
  }

  ret[["ll.lambda"]] <- ll.lambda

  # Convert list to data frame.
  empty <- sapply(ret, function(elem) length(elem) == 0)
  do.call(cbind.data.frame, ret[!empty])
}
