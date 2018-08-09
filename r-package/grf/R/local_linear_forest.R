#' Local Linear forest
#'
#' Trains a local linear forest that can be used to estimate
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
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed.
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
#' @return A trained local linear forest object.
#'
#' @examples \dontrun{
#' # Train a standard regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' forest = local_linear_forest(X, Y)
#' }
#'
#' @export
local_linear_forest <- function(X, Y,
                              sample.fraction = 0.5,
                              mtry = NULL, 
                              num.trees = 2000,
                              num.threads = NULL,
                              min.node.size = NULL,
                              honesty = TRUE,
                              alpha = NULL,
                              imbalance.penalty = NULL,
                              compute.oob.predictions = FALSE,
                              seed = NULL,
                              clusters = NULL,
                              samples_per_cluster = NULL,
                              tune.parameters = FALSE,
                              num.fit.trees = 10,
                              num.fit.reps = 100,
                              num.optimize.reps = 1000) {

    # return a standard regression forest with appropriate defaults
    forest = regression_forest(X, Y, sample.fraction,
                               mtry,
                               num.trees,
                               num.threads,
                               min.node.size,
                               honesty,
                               ci.group.size = 1,
                               alpha,
                               imbalance.penalty,
                               compute.oob.predictions,
                               seed,
                               clusters,
                               samples_per_cluster,
                               tune.parameters,
                               num.fit.trees,
                               num.fit.reps,
                               num.optimize.reps)

    forest
}

#' Predict with a local linear forest
#'
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param linear.correction.variables Optional subset of indexes for variables to be used in local
#'                   linear prediction. If left NULL, all variables are used.
#'                   We run a locally weighted linear regression on the included variables.
#'                   Please note that this is a beta feature still in development, and may slow down
#'                   prediction considerably. Defaults to NULL.
#' @param ll.lambda Ridge penalty for local linear predictions
#' @param tune.lambda Optional self-tuning for ridge penalty lambda. Defaults to FALSE.
#' @param lambda.path Optional list of lambdas to use for cross-validation, used if tune.lambda is TRUE.
#' @param ll.ridge.type Option to standardize ridge penalty by covariance ("standardized"),
#'                   or penalize all covariates equally ("identity").
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A vector of predictions.
#'
#' @examples \dontrun{
#' # Train the forest.
#' n = 50; p = 5
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' forest = local_linear_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' predictions = predict(forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' predictions.oob = predict(forest)
#' }
#'
#' @export
predict.local_linear_forest <- function(object, newdata = NULL,
                                      linear.correction.variables = NULL,
                                      ll.lambda = 0.01,
                                      tune.lambda = FALSE,
                                      lambda.path = NULL,
                                      ll.ridge.type = "standardized",
                                      num.threads = NULL,
                                      ...) {
    # force local linear prediction to be on
    X.orig = object[["X.orig"]]
    linear.correction.variables = validate_ll_vars(linear.correction.variables, ncol(X.orig))

    # then run regression forest predictions with local linear correction
    predictions = predict(object, newdata,
                          linear.correction.variables,
                          ll.lambda,
                          tune.lambda,
                          lambda.path,
                          ll.ridge.type,
                          nu.threads,
                          estimate.variance = FALSE)

    predictions
}
