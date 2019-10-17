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
#' @param sample.weights (experimental) Weights given to each observation in estimation.
#'                       If NULL, each observation receives equal weight. Default is NULL.
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
#' @param honesty.prune.leaves If true, prunes the estimation sample tree such that no leaves
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
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained instrumental forest object.
#' @export
instrumental_forest <- function(X, Y, W, Z,
                                Y.hat = NULL,
                                W.hat = NULL,
                                Z.hat = NULL,
                                num.trees = 2000,
                                sample.weights = NULL,
                                clusters = NULL,
                                samples.per.cluster = NULL,
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
                                compute.oob.predictions = TRUE,
                                num.threads = NULL,
                                seed = runif(1, 0, .Machine$integer.max)) {
  validate_X(X)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)
  Z <- validate_observations(Z, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_samples_per_cluster(samples.per.cluster, clusters)
  num.threads <- validate_num_threads(num.threads)

  if (!is.numeric(reduced.form.weight) | reduced.form.weight < 0 | reduced.form.weight > 1) {
    stop("Error: Invalid value for reduced.form.weight. Please give a value in [0,1].")
  }

  args.orthog = list(X = X,
                     num.trees = min(500, num.trees),
                     sample.weights = sample.weights,
                     clusters = clusters,
                     samples.per.cluster = samples.per.cluster,
                     sample.fraction = sample.fraction,
                     mtry = mtry,
                     min.node.size = 5,
                     honesty = TRUE,
                     honesty.fraction = 0.5,
                     honesty.prune.leaves = honesty.prune.leaves,
                     alpha = alpha,
                     imbalance.penalty = imbalance.penalty,
                     ci.group.size = 1,
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

  data <- create_data_matrices(X, outcome = Y - Y.hat, treatment = W - W.hat,
                               instrument = Z - Z.hat, sample.weights = sample.weights)
  forest <- instrumental_train(
    data$train.matrix, data$sparse.train.matrix,
    data$outcome.index, data$treatment.index, data$instrument.index, data$sample.weight.index,
    data$use.sample.weights,
    mtry, num.trees, min.node.size, sample.fraction, honesty, honesty.fraction,
    honesty.prune.leaves, ci.group.size, reduced.form.weight, alpha, imbalance.penalty, stabilize.splits, clusters,
    samples.per.cluster, compute.oob.predictions, num.threads, seed
  )

  class(forest) <- c("instrumental_forest", "grf")
  forest[["ci.group.size"]] <- ci.group.size
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["Z.orig"]] <- Z
  forest[["Y.hat"]] <- Y.hat
  forest[["W.hat"]] <- W.hat
  forest[["Z.hat"]] <- Z.hat
  forest[["clusters"]] <- clusters
  forest[["sample.weights"]] <- sample.weights
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
#' @method predict instrumental_forest
#' @export
predict.instrumental_forest <- function(object, newdata = NULL,
                                        num.threads = NULL,
                                        estimate.variance = FALSE,
                                        ...) {

  # If possible, use pre-computed predictions.
  if (is.null(newdata) & !estimate.variance & !is.null(object$predictions)) {
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

  train.data <- create_data_matrices(X, outcome = Y.centered, treatment = W.centered, instrument = Z.centered)

  if (!is.null(newdata)) {
    validate_newdata(newdata, object$X.orig)
    data <- create_data_matrices(newdata)
    ret <- instrumental_predict(
      forest.short, train.data$train.matrix, train.data$sparse.train.matrix,
      train.data$outcome.index, train.data$treatment.index, train.data$instrument.index,
      data$train.matrix, data$sparse.train.matrix, num.threads, estimate.variance
    )
  } else {
    ret <- instrumental_predict_oob(
      forest.short, train.data$train.matrix, train.data$sparse.train.matrix,
      train.data$outcome.index, train.data$treatment.index, train.data$instrument.index,
      num.threads, estimate.variance
    )
  }

  # Convert list to data frame.
  empty <- sapply(ret, function(elem) length(elem) == 0)
  do.call(cbind.data.frame, ret[!empty])
}
