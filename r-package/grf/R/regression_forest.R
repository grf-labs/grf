library(DiceKriging)

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
regression_forest <- function(X, Y, sample.fraction = 0.5, mtry = NULL,
                              num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                              honesty = TRUE, ci.group.size = 2, alpha = 0.05, imbalance.penalty = 0.0,
                              seed = NULL, clusters = NULL, samples_per_cluster = NULL) {
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }

    mtry <- validate_mtry(mtry, X)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)

    clusters <- validate_clusters(clusters, X)
    samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)

    data <- create_data_matrices(X, Y)
    outcome.index <- ncol(X) + 1

    forest <- regression_train(data$default, data$sparse, outcome.index, mtry,
        num.trees, num.threads, min.node.size, sample.fraction, seed, honesty,
        ci.group.size, alpha, imbalance.penalty, clusters, samples_per_cluster)

    forest[["ci.group.size"]] <- ci.group.size
    forest[["X.orig"]] <- X
    forest[["clusters"]] <- clusters

    class(forest) <- c("regression_forest", "grf")
    forest
}

#' Regression forest tuning
#'
#' Finds the optimal parameters to be used in training a regression forest. This method
#' currently tunes over min.node.size, sample.fraction, alpha, and imbalance.penalty.
#' Please see the method 'regression_forest' for a description of the standard forest
#' parameters. Note that if fixed values can be supplied for any of the parameters mentioned
#' above, and in that case, that parameter will not be tuned. For example, if this method is
#' called with min.node.size = 10 and alpha = 0.7, then those parameter values will be treated
#' as fixed, and only sample.fraction and imbalance.penalty will be tuned.
#'
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param num.fit.reps The number of forests used to fit the tuning model.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters.
#'
#' @return A list consisting of the optimal parameter values ('params') along with their debiased
#'         error ('error').
#'
#' @examples \dontrun{
#' # Find the optimal tuning parameters.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' params = tune_regression_forest(X, Y)$params
#'
#' # Use these parameters to train a regression forest.
#' tuned.forest = regression_forest(X, Y, num.trees = 1000,
#'     min.node.size = params["min.node.size"],
#'     sample.fraction = params["sample.fraction"],
#'     alpha = params["alpha"],
#'     imbalance.penalty = params["imbalance.penalty"])
#' }
#'
#' @export
tune_regression_forest <- function(X, Y,
                                   num.fit.trees = 10,
                                   num.fit.reps = 100,
                                   num.optimize.reps = 1000,
                                   min.node.size = NULL,
                                   sample.fraction = NULL,
                                   alpha = NULL,
                                   imbalance.penalty = NULL,
                                   num.threads = NULL,
                                   honesty = TRUE,
                                   seed = NULL,
                                   clusters = NULL,
                                   samples_per_cluster = NULL) {
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }

    num.threads <- validate_num_threads(num.threads)
    seed <- validate_seed(seed)
    clusters <- validate_clusters(clusters, X)
    samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)

    data <- create_data_matrices(X, Y)
    outcome.index <- ncol(X) + 1

    # Separate out the tuning parameters with supplied values, and those that were
    # left as 'NULL'. We will only tune those parameters that the user didn't supply.
    all.params = c(min.node.size = null_to_na(min.node.size),
      sample.fraction = null_to_na(sample.fraction),
      alpha = null_to_na(alpha),
      imbalance.penalty = null_to_na(imbalance.penalty))
    fixed.params = all.params[!is.na(all.params)]
    tuning.params = all.params[is.na(all.params)]

    # Train several mini-forests, and gather their debiased OOB error estimates.
    num.params = length(tuning.params)
    fit.draws = matrix(runif(num.fit.reps * num.params), num.fit.reps, num.params)
    colnames(fit.draws) = names(tuning.params)

    debiased.errors = apply(fit.draws, 1, function(draw) {
        params = c(fixed.params, get_params(X, draw))
        small.forest = regression_forest(X, Y,
                                         num.threads = num.threads, honesty = honesty, seed = seed,
                                         clusters = clusters, samples_per_cluster = samples_per_cluster,
                                         num.trees = num.fit.trees,
                                         min.node.size = params["min.node.size"],
                                         sample.fraction = params["sample.fraction"],
                                         alpha = params["alpha"],
                                         imbalance.penalty = params["imbalance.penalty"]
        )
        prediction = predict(small.forest)
        mean(prediction$debiased.error, na.rm = TRUE)
    })

    # Fit the 'dice kriging' model to these error estimates.
    # Note that in the 'km' call, the kriging package prints a large amount of information
    # about the fitting process. Here, capture its console output and discard it.
    variance.guess = rep(var(debiased.errors)/2, nrow(fit.draws))
    env = environment()
    capture.output(env$kriging.model <- km(design = data.frame(fit.draws),
                                           response = debiased.errors,
                                           noise.var = variance.guess))

    # To determine the optimal parameter values, predict using the kriging model at a large
    # number of random values, then select those that produced the lowest error.
    optimize.draws = matrix(runif(num.optimize.reps * num.params), num.optimize.reps, num.params)
    colnames(optimize.draws) = names(tuning.params)
    model.surface = predict(kriging.model, newdata=data.frame(optimize.draws), type = "SK")

    min.error = min(model.surface$mean)
    optimal.draw = optimize.draws[which.min(model.surface$mean),]
    optimal.params = get_params(X, optimal.draw)

    list("error"=min.error, "params"=optimal.params)
}

null_to_na <- function(param) {
  if (is.null(param)) NA else param
}

get_params <- function(X, draw) {
  result = c()
  for (param in names(draw)) {
    if (param == "min.node.size") {
      value = floor(2^(draw[param] * (log(nrow(X)) / log(2) - 4)))
    } else if (param == "sample.fraction") {
      value = 0.05 + 0.45 * draw[2]
    } else if (param == "alpha") {
      value = draw[3]/4
    } else if (param == "imbalance.penalty") {
      value = -log(draw[4])
    } else {
      stop("Unrecognized parameter name provided: ", param)
    }
    result = c(result, value)
  }
  result
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
#' @param local.linear Optional local linear prediction correction
#' @param lambda Ridge penalty for local linear predictions
#' @ridge.type Option to standardize ridge penalty by covariance ("standardized"),
#'                    or penalize all covariates equally ("identity").
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

    if(ridge.type == "standardized"){
        ridge_type = 0
    } else if(ridge.type == "identity"){
        ridge_type = 1
    } else{
        stop("Error: invalid ridge type")
    }

    forest.short <- object[-which(names(object) == "X.orig")]
    X.orig = object[["X.orig"]]

    if (!is.null(newdata)) {
        data <- create_data_matrices(newdata)

        if(!local.linear){
            regression_predict(forest.short, data$default, data$sparse,
                num.threads, ci.group.size)
        } else{
            training.data <- create_data_matrices(X.orig)
            local_linear_predict(forest.short, data$default, training.data$default, data$sparse,
                training.data$sparse, lambda, ridge_type, num.threads)
        }

    } else {
        data <- create_data_matrices(X.orig)

        if(!local.linear){
            regression_predict_oob(forest.short, data$default, data$sparse, num.threads, ci.group.size)
        } else{
            local_linear_predict_oob(forest.short, data$default, data$sparse, lambda, ridge_type,
                num.threads)
        }
    }
}
