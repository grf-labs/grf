#' Causal forest tuning
#'
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
#' @param Y.hat Estimates of the expected responses E[Y | Xi], marginalizing
#'              over treatment. See section 6.1.1 of the GRF paper for
#'              further discussion of this quantity.
#' @param W.hat Estimates of the treatment propensities E[W | Xi].
#' @param sample.weights (experimental) Weights given to an observation in estimation.
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
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#'                         determining the imbalance of a split. Default is TRUE.
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2.
#' @param tune.parameters A vector of parameter names to tune.
#'  If "all": all tunable parameters are tuned by cross-validation. The following parameters are
#'  tunable: ("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
#'   "honesty.prune.leaves", "alpha", "imbalance.penalty"). If honesty is FALSE the honesty.* parameters are not tuned.
#'  Default is "all".
#' @param tune.num.trees The number of trees in each 'mini forest' used to fit the tuning model. Default is 50.
#' @param tune.num.reps The number of forests used to fit the tuning model. Default is 100.
#' @param tune.num.draws The number of random parameter values considered when using the model
#'                          to select the optimal parameters. Default is 1000.
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
#' n <- 500
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
#' @export
tune_causal_forest <- function(X, Y, W, Y.hat, W.hat,
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
                              stabilize.splits = TRUE,
                              ci.group.size = 2,
                              tune.parameters = "all",
                              tune.num.trees = 200,
                              tune.num.reps = 50,
                              tune.num.draws = 1000,
                              num.threads = NULL,
                              seed = runif(1, 0, .Machine$integer.max)) {
  validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)

  all.tunable.params <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
                          "honesty.prune.leaves", "alpha", "imbalance.penalty")

  default.parameters <- list(sample.fraction = 0.5,
                             mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                             min.node.size = 5,
                             honesty.fraction = 0.5,
                             honesty.prune.leaves = TRUE,
                             alpha = 0.05,
                             imbalance.penalty = 0)

  Y.centered <- Y - Y.hat
  W.centered <- W - W.hat
  data <- create_data_matrices(X, outcome = Y.centered, treatment = W.centered,
                              sample.weights = sample.weights)
  nrow.X <- nrow(X)
  ncol.X <- ncol(X)
  args <- list(clusters = clusters,
               samples.per.cluster = samples.per.cluster,
               sample.fraction = sample.fraction,
               mtry = mtry,
               min.node.size = min.node.size,
               honesty = honesty,
               honesty.fraction = honesty.fraction,
               honesty.prune.leaves = honesty.prune.leaves,
               alpha = alpha,
               stabilize.splits = stabilize.splits,
               imbalance.penalty = imbalance.penalty,
               ci.group.size = ci.group.size,
               num.threads = num.threads,
               seed = seed,
               reduced.form.weight = 0)

  if (identical(tune.parameters, "all")) {
    tune.parameters <- all.tunable.params
  } else {
    tune.parameters <- unique(match.arg(tune.parameters, all.tunable.params, several.ok = TRUE))
  }
  if (!honesty) {
    tune.parameters <- tune.parameters[!grepl("honesty", tune.parameters)]
  }

  tune.parameters.defaults <- default.parameters[tune.parameters]
  train <- causal_train

  tuning.output <- tune_forest(data = data,
                               nrow.X = nrow.X,
                               ncol.X = ncol.X,
                               args = args,
                               tune.parameters = tune.parameters,
                               tune.parameters.defaults = tune.parameters.defaults,
                               num.fit.trees = tune.num.trees,
                               num.fit.reps = tune.num.reps,
                               num.optimize.reps = tune.num.draws,
                               train = train)

  tuning.output
}
