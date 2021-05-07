#' Boosted regression forest (experimental)
#'
#' Trains a boosted regression forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]. Selects
#' number of boosting iterations based on cross-validation. This functionality
#' is experimental and will likely change in future releases.
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param sample.weights Weights given to each observation in estimation.
#'                       If NULL, each observation receives the same weight. Default is NULL.
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
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2.
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if FALSE
#'                        NULL parameters are set to defaults. Default is FALSE.
#' @param tune.num.trees The number of trees in each 'mini forest' used to fit the tuning model. Default is 10.
#' @param tune.num.reps The number of forests used to fit the tuning model. Default is 100.
#' @param tune.num.draws The number of random parameter values considered when using the model
#'                          to select the optimal parameters. Default is 1000.
#' @param boost.steps The number of boosting iterations. If NULL, selected by cross-validation. Default is NULL.
#' @param boost.error.reduction If boost.steps is NULL, the percentage of previous steps' error that must be estimated
#'                  by cross validation in order to take a new step, default 0.97.
#' @param boost.max.steps The maximum number of boosting iterations to try when boost.steps=NULL. Default is 5.
#' @param boost.trees.tune If boost.steps is NULL, the number of trees used to test a new boosting step when tuning
#'        boost.steps. Default is 10.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param seed The seed for the C++ random number generator.
#'
#' @return A boosted regression forest object. $error contains the mean debiased error for each step, and $forests
#'         contains the trained regression forest for each step.
#'
#' @examples
#' \donttest{
#' # Train a boosted regression forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' boosted.forest <- boosted_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' boost.pred <- predict(boosted.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' boost.pred <- predict(boosted.forest)
#'
#' # Check how many boosting iterations were used
#' print(length(boosted.forest$forests))
#' }
#'
#' @export
boosted_regression_forest <- function(X, Y,
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
                                      ci.group.size = 2,
                                      tune.parameters = "none",
                                      tune.num.trees = 10,
                                      tune.num.reps = 100,
                                      tune.num.draws = 1000,
                                      boost.steps = NULL,
                                      boost.error.reduction = 0.97,
                                      boost.max.steps = 5,
                                      boost.trees.tune = 10,
                                      num.threads = NULL,
                                      seed = runif(1, 0, .Machine$integer.max)) {
  boost.error.reduction <- validate_boost_error_reduction(boost.error.reduction)
  boosted.forest <- NULL
  boosted.forest[["forests"]] <- list()
  boosted.forest[["error"]] <- list()
  forest.Y <- regression_forest(X, Y,
    sample.weights = sample.weights,
    sample.fraction = sample.fraction,
    mtry = mtry, tune.parameters = tune.parameters,
    tune.num.trees = tune.num.trees, tune.num.reps = tune.num.reps,
    tune.num.draws = tune.num.draws, num.trees = num.trees,
    num.threads = num.threads,
    min.node.size = min.node.size, honesty = honesty,
    honesty.fraction = honesty.fraction,
    honesty.prune.leaves = honesty.prune.leaves,
    seed = seed, ci.group.size = ci.group.size,
    alpha = alpha,
    imbalance.penalty = imbalance.penalty,
    clusters = clusters, equalize.cluster.weights = equalize.cluster.weights
  )
  current.pred <- predict(forest.Y, num.threads = num.threads)
  # save tuned parameters for use on future boosting iterations
  tunable.params <- forest.Y$tunable.params
  Y.hat <- current.pred$predictions
  error.debiased <- current.pred$debiased.error
  boosted.forest[["forests"]][[1]] <- forest.Y
  boosted.forest[["error"]][[1]] <- mean(error.debiased)

  step <- 1

  while (step <- step + 1) {
    Y.resid <- Y - Y.hat
    # do termination checks
    if (!is.null(boost.steps)) {
      if (step > boost.steps) {
        break
      }
    } else if (step > boost.max.steps) {
      break
    } else {
      # do cross validation check
      forest.small <- regression_forest(X, Y.resid,
        sample.weights = sample.weights,
        sample.fraction = as.numeric(tunable.params["sample.fraction"]),
        mtry = as.numeric(tunable.params["mtry"]), tune.parameters = "none",
        num.trees = boost.trees.tune,
        num.threads = num.threads,
        min.node.size = as.numeric(tunable.params["min.node.size"]),
        honesty = honesty,
        honesty.fraction = as.numeric(tunable.params["honesty.fraction"]),
        honesty.prune.leaves = as.numeric(tunable.params["honesty.prune.leaves"]),
        seed = seed, ci.group.size = ci.group.size,
        alpha = as.numeric(tunable.params["alpha"]),
        imbalance.penalty = as.numeric(tunable.params["imbalance.penalty"]),
        clusters = clusters, equalize.cluster.weights = equalize.cluster.weights
      )
      step.error.approx <- predict(forest.small, num.threads = num.threads)$debiased.error
      if (!(mean(step.error.approx, na.rm = TRUE) <= boost.error.reduction * mean(error.debiased, na.rm = TRUE))) {
        break
      }
    }

    forest.resid <- regression_forest(X, Y.resid,
      sample.weights = sample.weights,
      sample.fraction = as.numeric(tunable.params["sample.fraction"]),
      mtry = as.numeric(tunable.params["mtry"]), tune.parameters = "none",
      num.trees = num.trees,
      num.threads = num.threads,
      min.node.size = as.numeric(tunable.params["min.node.size"]),
      honesty = honesty,
      honesty.fraction = as.numeric(tunable.params["honesty.fraction"]),
      honesty.prune.leaves = as.numeric(tunable.params["honesty.prune.leaves"]),
      seed = seed, ci.group.size = ci.group.size,
      alpha = as.numeric(tunable.params["alpha"]),
      imbalance.penalty = as.numeric(tunable.params["imbalance.penalty"]),
      clusters = clusters, equalize.cluster.weights = equalize.cluster.weights
    )

    current.pred <- predict(forest.resid, num.threads = num.threads)
    Y.hat <- Y.hat + current.pred$predictions
    error.debiased <- current.pred$debiased.error
    boosted.forest[["forests"]][[step]] <- forest.resid
    boosted.forest[["error"]][[step]] <- mean(error.debiased)
  }
  boosted.forest[["predictions"]] <- Y.hat
  class(boosted.forest) <- c("boosted_regression_forest")
  boosted.forest
}

#' Predict with a boosted regression forest.
#'
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order
#' @param boost.predict.steps Number of boosting iterations to use for prediction. If blank, uses the full number of
#'        steps for the object given
#' @param num.threads the number of threads used in prediction
#' @param ... Additional arguments (currently ignored).
#'
#' @return A vector of predictions.
#'
#' @examples
#' \donttest{
#' # Train a boosted regression forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' r.boosted.forest <- boosted_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' r.pred <- predict(r.boosted.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred <- predict(r.boosted.forest)
#' }
#'
#' @method predict boosted_regression_forest
#' @export
predict.boosted_regression_forest <- function(object, newdata = NULL,
                                              boost.predict.steps = NULL,
                                              num.threads = NULL,
                                              ...) {

  # If not on new data, use pre-computed predictions
  if (is.null(newdata)) {
    return(data.frame(predictions = object$predictions))
  } else {
    forests <- object[["forests"]]
    if (is.null(boost.predict.steps)) {
      boost.predict.steps <- length(forests)
    } else {
      boost.predict.steps <- min(boost.predict.steps, length(forests))
    }
    Y.hat <- 0
    for (f in 1:boost.predict.steps) {
      Y.hat <- Y.hat + predict(forests[[f]], newdata, num.threads = num.threads)$predictions
    }
  }
  data.frame(predictions = Y.hat)
}
