#' Causal forest
#' 
#' Trains a causal forest that can be used to estimate
#' conditional average treatment effects tau(X). When
#' the treatment assignment W is binary and unconfounded,
#' we have tau(X) = E[Y(1) - Y(0) | X = x], where Y(0) and
#' Y(1) are potential outcomes corresponding to the two possible
#' treatment states. When W is continuous, we effectively estimate
#' an average partial effect Cov[Y, W | X = x] / Var[W | X = x],
#' and interpret it as a treatment effect given unconfoundedness.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
#' @param Y.hat Estimates of the expected responses E[Y | Xi], marginalizing
#'              over treatment. If Y.hat = NULL, these are estimated using
#'              a separate regression forest. See section 6.1.1 of the GRF paper for
#'              further discussion of this quantity.
#' @param W.hat Estimates of the treatment propensities E[W | Xi]. If W.hat = NULL,
#'              these are estimated using a separate regression forest.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting).
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds 
#'                         to set J1 in the notation of the paper. When using the defaults (honesty = TRUE and 
#'                         honesty.fraction = NULL), half of the data will be used for determining splits
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized.
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#'                         determining the imbalance of a split (experimental).
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed.
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
#'
#' # In some examples, pre-fitting models for Y and W separately may
#' # be helpful (e.g., if different models use different covariates).
#' # In some applications, one may even want to get Y.hat and W.hat
#' # using a completely different method (e.g., boosting).
#' n = 2000; p = 20
#' X = matrix(rnorm(n * p), n, p)
#' TAU = 1 / (1 + exp(-X[, 3]))
#' W = rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
#' Y = pmax(X[, 2] + X[, 3], 0) + rowMeans(X[, 4:6]) / 2 + W * TAU + rnorm(n)
#'
#' forest.W = regression_forest(X, W, tune.parameters = TRUE)
#' W.hat = predict(forest.W)$predictions
#'
#' forest.Y = regression_forest(X, Y, tune.parameters = TRUE)
#' Y.hat = predict(forest.Y)$predictions
#'
#' forest.Y.varimp = variable_importance(forest.Y)
#'
#' # Note: Forests may have a hard time when trained on very few variables
#' # (e.g., ncol(X) = 1, 2, or 3). We recommend not being too aggressive
#' # in selection.
#' selected.vars = which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)
#'
#' tau.forest = causal_forest(X[,selected.vars], Y, W,
#'                            W.hat = W.hat, Y.hat = Y.hat,
#'                            tune.parameters = TRUE)
#' tau.hat = predict(tau.forest)$predictions
#' }
#'
#' @export
causal_forest <- function(X, Y, W,
                          Y.hat = NULL,
                          W.hat = NULL,
                          sample.fraction = 0.5,
                          mtry = NULL,
                          num.trees = 2000,
                          num.threads = NULL,
                          min.node.size = NULL,
                          honesty = TRUE,
                          honesty.fraction = NULL,
                          ci.group.size = 2,
                          alpha = NULL,
                          imbalance.penalty = NULL,
                          stabilize.splits = TRUE,
                          compute.oob.predictions = TRUE,
                          seed = NULL,
                          clusters = NULL,
                          samples_per_cluster = NULL,
                          tune.parameters = FALSE,
                          num.fit.trees = 200,
                          num.fit.reps = 50,
                          num.optimize.reps = 1000) {
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    if(length(W) != nrow(X)) { stop("W has incorrect length.") }
    
    num.threads <- validate_num_threads(num.threads)
    seed <- validate_seed(seed)
    clusters <- validate_clusters(clusters, X)
    samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)
    honesty.fraction <- validate_honesty_fraction(honesty.fraction, honesty)
    
    reduced.form.weight <- 0

    if (is.null(Y.hat)) {
      forest.Y <- regression_forest(X, Y, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
      Y.hat <- predict(forest.Y)$predictions
    } else if (length(Y.hat) == 1) {
      Y.hat <- rep(Y.hat, nrow(X))
    } else if (length(Y.hat) != nrow(X)) {
      stop("Y.hat has incorrect length.")
    }

    if (is.null(W.hat)) {
      forest.W <- regression_forest(X, W, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
      W.hat <- predict(forest.W)$predictions
    } else if (length(W.hat) == 1) {
      W.hat <- rep(W.hat, nrow(X))
    } else if (length(W.hat) != nrow(X)) {
      stop("W.hat has incorrect length.")
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
                                          honesty.fraction = honesty.fraction,
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
                                 coerce_honesty_fraction(honesty.fraction),
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
    forest[["tunable.params"]] <- tunable.params

    class(forest) <- c("causal_forest", "grf")

    if (compute.oob.predictions) {
        oob.pred <- predict(forest)
        forest[["predictions"]] <- oob.pred$predictions
        forest[["debiased.error"]] <- oob.pred$debiased.error
    }

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
#' @method predict causal_forest
#' @export
predict.causal_forest <- function(object, newdata = NULL, num.threads = NULL, estimate.variance = FALSE, ...) {

    # If possible, use pre-computed predictions.
    if (is.null(newdata) & !estimate.variance & !is.null(object$predictions)) {
        return(data.frame(predictions=object$predictions,
                          debiased.error=object$debiased.error))
    }

    num.threads <- validate_num_threads(num.threads)

    if (estimate.variance) {
        ci.group.size = object$ci.group.size
    } else {
        ci.group.size = 1
    }

    forest.short <- object[-which(names(object) == "X.orig")]
    if (!is.null(newdata)) {
        data <- create_data_matrices(newdata)
        ret <- instrumental_predict(forest.short, data$default, data$sparse,
                                    num.threads, ci.group.size)
    } else {
        data <- create_data_matrices(object[["X.orig"]])
        ret <- instrumental_predict_oob(forest.short, data$default, data$sparse,
                                        num.threads, ci.group.size)
    }

    # Convert list to data frame.
    empty = sapply(ret, function(elem) length(elem) == 0)
    do.call(cbind.data.frame, ret[!empty])
}

