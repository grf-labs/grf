#' Causal forest
#' 
#' Trains a causal forest that can be used to estimate
#' conditional average treatment effects tau(X). When
#' the treatment assignmnet W is binary and unconfounded,
#' we have tau(X) = E[Y(1) - Y(0) | X = x], where Y(0) and
#' Y(1) are potential outcomes corresponding to the two possible
#' treatment states. When W is continuous, we effectively estimate
#' an average partical effect Cov[Y, W | X = x] / Var[W | X = x],
#' and interpret it as a treatment effect given unconfoundedness.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty is used, these subsamples will
#'                        further be cut in half.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#' @param honesty Whether or not honest splitting (i.e., sub-sample splitting) should be used.
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2.
#' @param precompute.nuisance Should we first run regression forests to estimate
#'                            y(x) = E[Y|X=x] and w(x) = E[W|X=x], and then run a
#'                            causal forest on the residuals? This approach is
#'                            recommended, computational resources permitting.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized.
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#'                         determining the imbalance of a split (experimental).
#' @param seed The seed of the C++ random number generator.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#' @param samples_per_cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster. Must be less than the size of the smallest cluster. If set to NULL
#'                            software will set this value to the size of the smallest cluster.#'
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if false
#'                        NULL parameters are set to defaults.
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param num.fit.reps The number of forests used to fit the tuning model.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters.
#'
#' @return A trained causal forest object.
#'
#' @examples \dontrun{
#' # Train a causal forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' c.forest = causal_forest(X, Y, W)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' c.pred = predict(c.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' c.forest = causal_forest(X, Y, W, num.trees = 4000)
#' c.pred = predict(c.forest, X.test, estimate.variance = TRUE)
#' }
#'
#' @export
causal_forest <- function(X, Y, W,
                          sample.fraction = 0.5,
                          mtry = NULL,
                          num.trees = 2000,
                          num.threads = NULL,
                          min.node.size = NULL,
                          honesty = TRUE,
                          ci.group.size = 2,
                          precompute.nuisance = TRUE,
                          alpha = NULL,
                          imbalance.penalty = NULL,
                          stabilize.splits = TRUE,
                          seed = NULL,
                          clusters = NULL,
                          samples_per_cluster = NULL,
                          tune.parameters = FALSE,
                          num.fit.trees = 40,
                          num.fit.reps = 100,
                          num.optimize.reps = 1000) {
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    if(length(W) != nrow(X)) { stop("W has incorrect length.") }
    
    num.threads <- validate_num_threads(num.threads)
    seed <- validate_seed(seed)
    clusters <- validate_clusters(clusters, X)
    samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)
    
    reduced.form.weight <- 0

    if (!precompute.nuisance) {
      Y.hat = 0
      W.hat = 0
    } else {
      forest.Y <- regression_forest(X, Y, sample.fraction = sample.fraction, mtry = mtry, 
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, 
                                    honesty = TRUE, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
      Y.hat <- predict(forest.Y)$predictions

      forest.W <- regression_forest(X, W, sample.fraction = sample.fraction, mtry = mtry, 
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, 
                                    honesty = TRUE, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
      W.hat <- predict(forest.W)$predictions
    }

    Y.centered = Y - Y.hat
    W.centered = W - W.hat

    if (tune.parameters) {
      tuning.output <- tune_causal_forest(X, Y.centered, W.centered,
                                          num.fit.trees = num.fit.trees,
                                          num.fit.reps = num.fit.reps,
                                          num.optimize.reps = num.optimize.reps,
                                          min.node.size = min.node.size,
                                          sample.fraction = sample.fraction,
                                          mtry = mtry,
                                          alpha = alpha,
                                          imbalance.penalty = imbalance.penalty,
                                          stabilize.splits = stabilize.splits,
                                          num.threads = num.threads,
                                          honesty = honesty,
                                          seed = seed,
                                          clusters = clusters,
                                          samples_per_cluster = samples_per_cluster)
      tunable.params <- tuning.output$params
    } else {
      tunable.params <- c(
        min.node.size = validate_min_node_size(min.node.size),
        sample.fraction = validate_sample_fraction(sample.fraction),
        mtry = validate_mtry(mtry, X),
        alpha = validate_alpha(alpha),
        imbalance.penalty = validate_imbalance_penalty(imbalance.penalty))
    }

    data <- create_data_matrices(X, Y.centered, W.centered)
    outcome.index <- ncol(X) + 1
    treatment.index <- ncol(X) + 2
    instrument.index <- treatment.index

    forest <- instrumental_train(data$default, data$sparse,
                                 outcome.index, treatment.index, instrument.index,
                                 as.numeric(tunable.params["mtry"]),
                                 num.trees,
                                 num.threads,
                                 as.numeric(tunable.params["min.node.size"]),
                                 as.numeric(tunable.params["sample.fraction"]),
                                 seed,
                                 honesty,
                                 ci.group.size,
                                 reduced.form.weight,
                                 as.numeric(tunable.params["alpha"]),
                                 as.numeric(tunable.params["imbalance.penalty"]),
                                 stabilize.splits,
                                 clusters,
                                 samples_per_cluster)

    forest[["ci.group.size"]] <- ci.group.size
    forest[["X.orig"]] <- X
    forest[["Y.orig"]] <- Y
    forest[["W.orig"]] <- W
    forest[["Y.hat"]] <- Y.hat
    forest[["W.hat"]] <- W.hat
    forest[["clusters"]] <- clusters
    
    class(forest) <- c("causal_forest", "grf")
    forest
}

#' Predict with a causal forest
#' 
#' Gets estimates of tau(x) using a trained causal forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions, along with (optional) variance estimates.
#'
#' @examples \dontrun{
#' # Train a causal forest.
#' n = 100; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' c.forest = causal_forest(X, Y, W)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' c.pred = predict(c.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' c.forest = causal_forest(X, Y, W, num.trees = 500)
#' c.pred = predict(c.forest, X.test, estimate.variance = TRUE)
#' }
#'
#' @export
predict.causal_forest <- function(object, newdata = NULL, num.threads = NULL, estimate.variance = FALSE, ...) {
  num.threads <- validate_num_threads(num.threads) 
  if (estimate.variance) {
      ci.group.size = object$ci.group.size
  } else {
      ci.group.size = 1
  }

  forest.short <- object[-which(names(object) == "X.orig")]
  if (!is.null(newdata)) {
      data <- create_data_matrices(newdata)
      instrumental_predict(forest.short, data$default, data$sparse,
                           num.threads, ci.group.size)
  } else {
      data <- create_data_matrices(object[["X.orig"]])
      instrumental_predict_oob(forest.short, data$default, data$sparse,
                               num.threads, ci.group.size)
  }
}

