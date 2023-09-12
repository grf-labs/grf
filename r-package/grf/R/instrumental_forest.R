#' Intrumental forest
#'
#' Trains an instrumental forest that can be used to estimate
#' conditional local average treatment effects tau(X) identified
#' using instruments. Formally, the forest estimates
#' tau(X) = Cov[Y, Z | X = x] / Cov[W, Z | X = x].
#' Note that when the instrument Z and treatment assignment W
#' coincide, an instrumental forest is equivalent to a causal forest.
#'
#' @param X The covariates used in the instrumental regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
#' @param Z The instrument (may be binary or real).
#' @param Y.hat Estimates of the expected responses E[Y | Xi], marginalizing
#'              over treatment. If Y.hat = NULL, these are estimated using
#'              a separate regression forest. Default is NULL.
#' @param W.hat Estimates of the treatment propensities E[W | Xi]. If W.hat = NULL,
#'              these are estimated using a separate regression forest. Default is NULL.
#' @param Z.hat Estimates of the instrument propensities E[Z | Xi]. If Z.hat = NULL,
#'              these are estimated using a separate regression forest. Default is NULL.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param sample.weights Weights given to each observation in estimation.
#'                       If NULL, each observation receives equal weight. Default is NULL.
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
#' @param stabilize.splits Whether or not the instrument should be taken into account when
#'                         determining the imbalance of a split. Default is TRUE.
#' @param ci.group.size The forst will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2.
#' @param reduced.form.weight Whether splits should be regularized towards a naive
#'                            splitting criterion that ignores the instrument (and
#'                            instead emulates a causal forest).
#' @param tune.parameters A vector of parameter names to tune.
#'  If "all": all tunable parameters are tuned by cross-validation. The following parameters are
#'  tunable: ("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
#'   "honesty.prune.leaves", "alpha", "imbalance.penalty"). If honesty is FALSE the honesty.* parameters are not tuned.
#'  Default is "none" (no parameters are tuned).
#' @param tune.num.trees The number of trees in each 'mini forest' used to fit the tuning model. Default is 200.
#' @param tune.num.reps The number of forests used to fit the tuning model. Default is 50.
#' @param tune.num.draws The number of random parameter values considered when using the model
#'                          to select the optimal parameters. Default is 1000.
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained instrumental forest object.
#'
#' @references Athey, Susan, Julie Tibshirani, and Stefan Wager. "Generalized Random Forests".
#'  Annals of Statistics, 47(2), 2019.
#'
#' @examples
#' \donttest{
#' # Train an instrumental forest.
#' n <- 2000
#' p <- 5
#' X <- matrix(rbinom(n * p, 1, 0.5), n, p)
#' Z <- rbinom(n, 1, 0.5)
#' Q <- rbinom(n, 1, 0.5)
#' W <- Q * Z
#' tau <-  X[, 1] / 2
#' Y <- rowSums(X[, 1:3]) + tau * W + Q + rnorm(n)
#' iv.forest <- instrumental_forest(X, Y, W, Z)
#'
#' # Predict on out-of-bag training samples.
#' iv.pred <- predict(iv.forest)
#'
#' # Estimate a (local) average treatment effect.
#' average_treatment_effect(iv.forest)
#' }
#'
#' @export
instrumental_forest <- function(X, Y, W, Z,
                                Y.hat = NULL,
                                W.hat = NULL,
                                Z.hat = NULL,
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
                                stabilize.splits = TRUE,
                                ci.group.size = 2,
                                reduced.form.weight = 0,
                                tune.parameters = "none",
                                tune.num.trees = 200,
                                tune.num.reps = 50,
                                tune.num.draws = 1000,
                                compute.oob.predictions = TRUE,
                                num.threads = NULL,
                                seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)
  Z <- validate_observations(Z, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)

  if (!is.numeric(reduced.form.weight) | reduced.form.weight < 0 | reduced.form.weight > 1) {
    stop("Error: Invalid value for reduced.form.weight. Please give a value in [0,1].")
  }

  all.tunable.params <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
                          "honesty.prune.leaves", "alpha", "imbalance.penalty")
  default.parameters <- list(sample.fraction = 0.5,
                             mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                             min.node.size = 5,
                             honesty.fraction = 0.5,
                             honesty.prune.leaves = TRUE,
                             alpha = 0.05,
                             imbalance.penalty = 0)

  args.orthog = list(X = X,
                     num.trees = min(500, num.trees),
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
                     ci.group.size = 1,
                     tune.parameters = tune.parameters,
                     num.threads = num.threads,
                     seed = seed)

  if (is.null(Y.hat)) {
    forest.Y <- do.call(regression_forest, c(Y = list(Y), args.orthog))
    Y.hat <- predict(forest.Y)$predictions
  } else if (length(Y.hat) == 1) {
    Y.hat <- rep(Y.hat, nrow(X))
  } else if (length(Y.hat) != nrow(X)) {
    stop("Y.hat has incorrect length.")
  }

  if (is.null(W.hat)) {
    forest.W <- do.call(regression_forest, c(Y = list(W), args.orthog))
    W.hat <- predict(forest.W)$predictions
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, nrow(X))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }

  if (is.null(Z.hat)) {
    forest.Z <- do.call(regression_forest, c(Y = list(Z), args.orthog))
    Z.hat <- predict(forest.Z)$predictions
  } else if (length(Z.hat) == 1) {
    Z.hat <- rep(Z.hat, nrow(X))
  } else if (length(Z.hat) != nrow(X)) {
    stop("Z.hat has incorrect length.")
  }

  data <- create_train_matrices(X, outcome = Y - Y.hat, treatment = W - W.hat,
                                instrument = Z - Z.hat, sample.weights = sample.weights)
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
              reduced.form.weight = reduced.form.weight,
              compute.oob.predictions = compute.oob.predictions,
              num.threads = num.threads,
              seed = seed)

  tuning.output <- NULL
  if (!identical(tune.parameters, "none")) {
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
                                 train = instrumental_train)

    args <- utils::modifyList(args, as.list(tuning.output[["params"]]))
  }

  forest <- do.call.rcpp(instrumental_train, c(data, args))
  class(forest) <- c("instrumental_forest", "grf")
  forest[["seed"]] <- seed
  forest[["ci.group.size"]] <- ci.group.size
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["Z.orig"]] <- Z
  forest[["Y.hat"]] <- Y.hat
  forest[["W.hat"]] <- W.hat
  forest[["Z.hat"]] <- Z.hat
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["sample.weights"]] <- sample.weights
  forest[["tunable.params"]] <- args[all.tunable.params]
  forest[["tuning.output"]] <- tuning.output
  forest[["has.missing.values"]] <- has.missing.values

  forest
}

#' Predict with an instrumental forest
#'
#' Gets estimates of tau(x) using a trained instrumental forest.
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
#' @return Vector of predictions, along with (optional) variance estimates.
#'
#' @examples
#' \donttest{
#' # Train an instrumental forest.
#' n <- 2000
#' p <- 5
#' X <- matrix(rbinom(n * p, 1, 0.5), n, p)
#' Z <- rbinom(n, 1, 0.5)
#' Q <- rbinom(n, 1, 0.5)
#' W <- Q * Z
#' tau <-  X[, 1] / 2
#' Y <- rowSums(X[, 1:3]) + tau * W + Q + rnorm(n)
#' iv.forest <- instrumental_forest(X, Y, W, Z)
#'
#' # Predict on out-of-bag training samples.
#' iv.pred <- predict(iv.forest)
#'
#' # Estimate a (local) average treatment effect.
#' average_treatment_effect(iv.forest)
#' }
#'
#' @method predict instrumental_forest
#' @export
predict.instrumental_forest <- function(object, newdata = NULL,
                                        num.threads = NULL,
                                        estimate.variance = FALSE,
                                        ...) {

  # If possible, use pre-computed predictions.
  if (is.null(newdata) && !estimate.variance && !is.null(object$predictions)) {
    return(data.frame(
      predictions = object$predictions,
      debiased.error = object$debiased.error
    ))
  }

  num.threads <- validate_num_threads(num.threads)
  forest.short <- object[-which(names(object) == "X.orig")]

  X <- object[["X.orig"]]
  Y.centered <- object[["Y.orig"]] - object[["Y.hat"]]
  W.centered <- object[["W.orig"]] - object[["W.hat"]]
  Z.centered <- object[["Z.orig"]] - object[["Z.hat"]]

  train.data <- create_train_matrices(X, outcome = Y.centered, treatment = W.centered, instrument = Z.centered)
  args <- list(forest.object = forest.short,
               num.threads = num.threads,
               estimate.variance = estimate.variance)

  if (!is.null(newdata)) {
    validate_newdata(newdata, X, allow.na = TRUE)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(instrumental_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(instrumental_predict_oob, c(train.data, args))
  }

  # Convert list to data frame.
  empty <- sapply(ret, function(elem) length(elem) == 0)
  do.call(cbind.data.frame, ret[!empty])
}
