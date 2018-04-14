#' Regression forest
#'
#' Trains a regression forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
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
#' @param alpha A tuning parameter that controls the maximum imbalance of a split.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized.
#' @param seed The seed for the C++ random number generator.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#' @param samples_per_cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster. Must be less than the size of the smallest cluster. If set to NULL
#'                            software will set this value to the size of the smallest cluster.
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if false
#'                        NULL parameters are set to defaults.
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param num.fit.reps The number of forests used to fit the tuning model.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters.
#'
#' @return A trained regression forest object.
#'
#' @examples \dontrun{
#' # Train a standard regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.forest = regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' r.pred = predict(r.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred = predict(r.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' r.forest = regression_forest(X, Y, num.trees = 100)
#' r.pred = predict(r.forest, X.test, estimate.variance = TRUE)
#' }
#'
#' @export
regression_forest <- function(X, Y,
                              sample.fraction = 0.5,
                              mtry = NULL, 
                              num.trees = 2000,
                              num.threads = NULL,
                              min.node.size = NULL,
                              honesty = TRUE,
                              ci.group.size = 2,
                              alpha = NULL,
                              imbalance.penalty = NULL,
                              seed = NULL,
                              clusters = NULL,
                              samples_per_cluster = NULL,
                              tune.parameters = FALSE,
                              num.fit.trees = 10,
                              num.fit.reps = 100,
                              num.optimize.reps = 1000) {
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }

    num.threads <- validate_num_threads(num.threads)
    seed <- validate_seed(seed)
    clusters <- validate_clusters(clusters, X)
    samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)
    
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
    
    data <- create_data_matrices(X, Y)
    outcome.index <- ncol(X) + 1

    forest <- regression_train(data$default, data$sparse, outcome.index,
                               as.numeric(tunable.params["mtry"]),
                               num.trees,
                               num.threads,
                               as.numeric(tunable.params["min.node.size"]),
                               as.numeric(tunable.params["sample.fraction"]),
                               seed,
                               honesty,
                               ci.group.size,
                               as.numeric(tunable.params["alpha"]),
                               as.numeric(tunable.params["imbalance.penalty"]),
                               clusters,
                               samples_per_cluster)
    
    forest[["ci.group.size"]] <- ci.group.size
    forest[["X.orig"]] <- X
    forest[["clusters"]] <- clusters
    forest[["tunable.params"]] <- tunable.params
    
    class(forest) <- c("regression_forest", "grf")
    forest
}

#' Predict with a regression forest
#'
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param local.linear Optional local linear prediction correction. If TRUE,
#'                code will run a locally weighted ridge regression at each test point.
#'                Note that this is a beta feature still in development, and may slow down
#'                prediction considerably.
#' @param lambda Ridge penalty for local linear predictions
#' @param ridge.type Option to standardize ridge penalty by covariance ("standardized"),
#'                   or penalize all covariates equally ("identity").
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A vector of predictions.
#'
#' @examples \dontrun{
#' # Train a standard regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.forest = regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' r.pred = predict(r.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred = predict(r.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' r.forest = regression_forest(X, Y, num.trees = 100)
#' r.pred = predict(r.forest, X.test, estimate.variance = TRUE)
#' }
#'
#' @export
predict.regression_forest <- function(object, newdata = NULL,
                                      local.linear = FALSE,
                                      lambda = 0.0,
                                      ridge.type = "standardized",
                                      num.threads = NULL,
                                      estimate.variance = FALSE,
                                      ...) {
    num.threads <- validate_num_threads(num.threads)

    if (estimate.variance) {
        ci.group.size = object$ci.group.size
    } else {
        ci.group.size = 1
    }

    if (ridge.type == "standardized") {
        use_unweighted_penalty = 0
    } else if (ridge.type == "identity") {
        use_unweighted_penalty = 1
    } else {
        stop("Error: invalid ridge type")
    }

    forest.short <- object[-which(names(object) == "X.orig")]
    X.orig = object[["X.orig"]]

    if (!is.null(newdata)) {
        data <- create_data_matrices(newdata)
        if (!local.linear) {
            regression_predict(forest.short, data$default, data$sparse,
                num.threads, ci.group.size)
        } else {
            training.data <- create_data_matrices(X.orig)
            local_linear_predict(forest.short, data$default, training.data$default, data$sparse,
                training.data$sparse, lambda, use_unweighted_penalty, num.threads)
        }
    } else {
        data <- create_data_matrices(X.orig)
        if (!local.linear) {
            regression_predict_oob(forest.short, data$default, data$sparse, num.threads, ci.group.size)
        } else {
            local_linear_predict_oob(forest.short, data$default, data$sparse, lambda, use_unweighted_penalty,
                num.threads)
        }
    }
}
