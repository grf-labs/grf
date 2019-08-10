#' Causal forest tuning
#'
#' Finds the optimal parameters to be used in training a regression forest. This method
#' currently tunes over min.node.size, mtry, sample.fraction, alpha, and imbalance.penalty.
#' Please see the method 'causal_forest' for a description of the standard causal forest
#' parameters. Note that if fixed values can be supplied for any of the parameters mentioned
#' above, and in that case, that parameter will not be tuned. For example, if this method is
#' called with min.node.size = 10 and alpha = 0.7, then those parameter values will be treated
#' as fixed, and only sample.fraction and imbalance.penalty will be tuned.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
#' @param Y.hat Estimates of the expected responses E[Y | Xi], marginalizing
#'              over treatment. See section 6.1.1 of the GRF paper for
#'              further discussion of this quantity.
#' @param W.hat Estimates of the treatment propensities E[W | Xi].
#' @param sample.weights Weights defining the population on which we want our estimator of tau(x) to perform well
#'                       on average. If NULL, this is the population from which X1 ... Xn are sampled. Otherwise,
#'                       it is a reweighted version, in which we observe Xi with probability proportional to
#'                       sample.weights[i]. Default is NULL.
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model. Default is 200.
#' @param num.fit.reps The number of forests used to fit the tuning model. Default is 50.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters. Default is 1000.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction. Default is 0.5.
#' @param mtry Number of variables tried for each split. Default is
#'             \eqn{\sqrt p + 20} where p is the number of variables.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#'                      Default is 5.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting). Default is TRUE.
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. When using the defaults (honesty = TRUE and
#'                         honesty.fraction = NULL), half of the data will be used for determining splits.
#'                         Default is 0.5.
#' @param prune.empty.leaves (experimental) If true, prunes the estimation sample tree such that no leaves
#'  are empty. If false, keep the same tree as determined in the splits sample (if an empty leave is encountered, that
#'  tree is skipped and does not contribute to the estimate). Setting this to false may improve performance on
#'  small/marginally powered data, but requires more trees. Only applies if honesty is enabled. Default is TRUE.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. Default is 0.05.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized. Default is 0.
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#'                         determining the imbalance of a split (experimental). Default is TRUE.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'                 Default is NULL (ignored).
#' @param samples.per.cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster. Must be less than the size of the smallest cluster. If set to NULL
#'                            software will set this value to the size of the smallest cluster. Default is NULL.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A list consisting of the optimal parameter values ('params') along with their debiased
#'         error ('error').
#'
#' @examples
#' \dontrun{
#' # Find the optimal tuning parameters.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' Y.hat <- predict(regression_forest(X, Y))$predictions
#' W.hat <- rep(0.5, n)
#' params <- tune_causal_forest(X, Y, W, Y.hat, W.hat)$params
#'
#' # Use these parameters to train a regression forest.
#' tuned.forest <- causal_forest(X, Y, W,
#'   Y.hat = Y.hat, W.hat = W.hat, num.trees = 1000,
#'   min.node.size = as.numeric(params["min.node.size"]),
#'   sample.fraction = as.numeric(params["sample.fraction"]),
#'   mtry = as.numeric(params["mtry"]),
#'   alpha = as.numeric(params["alpha"]),
#'   imbalance.penalty = as.numeric(params["imbalance.penalty"])
#' )
#' }
#'
#' @importFrom stats runif
#' @importFrom utils capture.output
#' @export
tune_causal_forest <- function(X, Y, W, Y.hat, W.hat,
                               sample.weights = NULL,
                               num.fit.trees = 100,
                               num.fit.reps = 50,
                               num.optimize.reps = 1000,
                               min.node.size = NULL,
                               sample.fraction = 0.5,
                               mtry = NULL,
                               alpha = NULL,
                               imbalance.penalty = NULL,
                               stabilize.splits = TRUE,
                               honesty = TRUE,
                               honesty.fraction = NULL,
                               prune.empty.leaves = TRUE,
                               clusters = NULL,
                               samples.per.cluster = NULL,
                               num.threads = NULL,
                               seed = NULL) {
  validate_X(X)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)

  num.threads <- validate_num_threads(num.threads)
  seed <- validate_seed(seed)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_samples_per_cluster(samples.per.cluster, clusters)
  ci.group.size <- 1
  reduced.form.weight <- 0
  honesty.fraction <- validate_honesty_fraction(honesty.fraction, honesty)

  data <- create_data_matrices(X, Y - Y.hat, W - W.hat, sample.weights = sample.weights)
  outcome.index <- ncol(X) + 1
  treatment.index <- ncol(X) + 2
  sample.weight.index <- ncol(X) + 3

  # Separate out the tuning parameters with supplied values, and those that were
  # left as 'NULL'. We will only tune those parameters that the user didn't supply.
  all.params <- get_initial_params(min.node.size, sample.fraction, mtry, alpha, imbalance.penalty)
  fixed.params <- all.params[!is.na(all.params)]
  tuning.params <- all.params[is.na(all.params)]

  if (length(tuning.params) == 0) {
    return(list("error" = NA, "params" = c(all.params)))
  }

  # Train several mini-forests, and gather their debiased OOB error estimates.
  num.params <- length(tuning.params)
  fit.draws <- matrix(runif(num.fit.reps * num.params), num.fit.reps, num.params)
  colnames(fit.draws) <- names(tuning.params)
  compute.oob.predictions <- TRUE

  small.forest.errors <- apply(fit.draws, 1, function(draw) {
    params <- c(fixed.params, get_params_from_draw(X, draw))
    small.forest <- causal_train(
      data$default, data$sparse,
      outcome.index, treatment.index, sample.weight.index,
      !is.null(sample.weights),
      as.numeric(params["mtry"]),
      num.fit.trees,
      as.numeric(params["min.node.size"]),
      as.numeric(params["sample.fraction"]),
      honesty,
      coerce_honesty_fraction(honesty.fraction),
      prune.empty.leaves,
      ci.group.size,
      reduced.form.weight,
      as.numeric(params["alpha"]),
      as.numeric(params["imbalance.penalty"]),
      stabilize.splits,
      clusters,
      samples.per.cluster,
      compute.oob.predictions,
      num.threads,
      seed
    )
    prediction <- causal_predict_oob(
      small.forest, data$default, data$sparse,
      outcome.index, treatment.index, num.threads, FALSE
    )
    mean(prediction$debiased.error, na.rm = TRUE)
  })

  if (all(is.na(small.forest.errors))) {
    warning(paste0("Could not tune causal forest because all small forest error estimates were NA.\n",
                   "Consider increasing argument num.fit.trees."))
    out <- get_tuning_output(params = c(all.params), status = "failure")
    return(out)
  }

  if (sd(small.forest.errors) < 1e-10) {
    warning(paste0("Could not tune causal forest because small forest errors were nearly constant.\n",
                   "Consider increasing argument num.fit.trees."))
    out <- get_tuning_output(params = c(all.params), status = "failure")
    return(out)
  }

  # Fit the 'dice kriging' model to these error estimates.
  # Note that in the 'km' call, the kriging package prints a large amount of information
  # about the fitting process. Here, capture its console output and discard it.
  variance.guess <- rep(var(small.forest.errors) / 2, nrow(fit.draws))
  kriging.model <- tryCatch({
    capture.output(
      DiceKriging::km(
        design = data.frame(fit.draws),
        response = small.forest.errors,
        noise.var = variance.guess
      ))
    },
    error = function(e) {
      warning(paste0("Dicekriging threw the following error: \n", e)
      NA
  })
  if (is.na(kriging.model)) {
    warning("Tuning was attempted but failed. Reverting to default parameters.")
    out <- get_tuning_output(params = pre.tuning.parameters, status = "failure")
    return(out)
  }

  # To determine the optimal parameter values, predict using the kriging model at a large
  # number of random values, then select those that produced the lowest error.
  optimize.draws <- matrix(runif(num.optimize.reps * num.params), num.optimize.reps, num.params)
  colnames(optimize.draws) <- names(tuning.params)
  model.surface <- predict(kriging.model, newdata = data.frame(optimize.draws), type = "SK")$mean
  tuned.params <- get_params_from_draw(X, optimize.draws)

  grid <- cbind(error = c(model.surface), tuned.params)
  small.forest.optimal.draw <- which.min(grid[, "error"])
  small.forest.optimal.params <- grid[small.forest.optimal.draw, -1]

  # To avoid the possibility of selection bias, re-train a moderately-sized forest
  # at the value chosen by the method above
  retrained.forest.params <- c(fixed.params, grid[small.forest.optimal.draw, -1])
  retrained.forest.num.trees <- num.fit.trees * 4
  retrained.forest <- causal_train(data$default, data$sparse,
                                  outcome.index, treatment.index, sample.weight.index,
                                  !is.null(sample.weights),
                                  retrained.forest.params["mtry"],
                                  retrained.forest.num.trees,
                                  retrained.forest.params["min.node.size"],
                                  retrained.forest.params["sample.fraction"],
                                  honesty,
                                  coerce_honesty_fraction(honesty.fraction),
                                  prune.empty.leaves,
                                  ci.group.size,
                                  reduced.form.weight,
                                  retrained.forest.params["alpha"],
                                  retrained.forest.params["imbalance.penalty"],
                                  stabilize.splits,
                                  clusters,
                                  samples.per.cluster,
                                  compute.oob.predictions,
                                  num.threads,
                                  seed)

  retrained.forest.prediction <- causal_predict_oob(
    retrained.forest, data$default, data$sparse,
    outcome.index, treatment.index, num.threads, FALSE)

  retrained.forest.error <- mean(retrained.forest.prediction$debiased.error, na.rm = TRUE)


  # Train a forest with default parameters, and check its predicted error.
  # This improves our chances of not doing worse than default
  default.params <- c(
    min.node.size = validate_min_node_size(min.node.size),
    sample.fraction = validate_sample_fraction(sample.fraction),
    mtry = validate_mtry(mtry, X),
    alpha = validate_alpha(alpha),
    imbalance.penalty = validate_imbalance_penalty(imbalance.penalty)
  )
  default.params[!is.na(all.params)] <- all.params[!is.na(all.params)]
  num.default.forest.trees <- num.fit.trees * 4

  default.forest <- causal_train(
    data$default, data$sparse,
    outcome.index, treatment.index, sample.weight.index,
    !is.null(sample.weights),
    default.params["mtry"],
    num.default.forest.trees,
    default.params["min.node.size"],
    default.params["sample.fraction"],
    honesty,
    coerce_honesty_fraction(honesty.fraction),
    prune.empty.leaves,
    ci.group.size,
    reduced.form.weight,
    default.params["alpha"],
    default.params["imbalance.penalty"],
    stabilize.splits,
    clusters,
    samples.per.cluster,
    compute.oob.predictions,
    num.threads,
    seed)

  default.forest.prediction <- causal_predict_oob(
    default.forest, data$default, data$sparse,
    outcome.index, treatment.index, num.threads, FALSE)

  default.forest.error <- mean(default.forest.prediction$debiased.error, na.rm = TRUE)

  # Now compare predicted errors: default parameters vs retrained parameter
  if (default.forest.error < retrained.forest.error) {
    out <- get_tuning_output(error = default.forest.error,
                             params = default.params,
                             grid = NA,
                             status = "default")
  } else {
    out <- get_tuning_output(error = retrained.forest.error,
                             params = retrained.forest.params,
                             grid = grid,
                             status = "tuned")
  }

  out
}
