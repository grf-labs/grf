#' Local linear forest tuning
#'
#' Finds the optimal ridge penalty for local linear causal prediction.
#'
#' @param forest The forest used for prediction.
#' @param linear.correction.variables Variables to use for local linear prediction. If left null,
#'          all variables are used. Default is NULL.
#' @param ll.weight.penalty Option to standardize ridge penalty by covariance (TRUE),
#'                            or penalize all covariates equally (FALSE). Defaults to FALSE.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param lambda.path Optional list of lambdas to use for cross-validation.
#' @return A list of lambdas tried, corresponding errors, and optimal ridge penalty lambda.
#'
#' @examples
#' \donttest{
#' # Find the optimal tuning parameters.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#'
#' forest <- causal_forest(X, Y, W)
#' tuned.lambda <- tune_ll_causal_forest(forest)
#'
#' # Use this parameter to predict from a local linear causal forest.
#' predictions <- predict(forest, linear.correction.variables = 1:p,
#'                        ll.lambda = tuned.lambda$lambda.min)
#' }
#'
#' @export
tune_ll_causal_forest <- function(forest,
                                  linear.correction.variables = NULL,
                                  ll.weight.penalty = FALSE,
                                  num.threads = NULL,
                                  lambda.path = NULL) {
  forest.short <- forest[-which(names(forest) == "X.orig")]
  X <- forest[["X.orig"]]
  Y <- forest[["Y.orig"]]
  W <- forest[["W.orig"]]
  Y.hat <- forest[["Y.hat"]]
  W.hat <- forest[["W.hat"]]
  Y.centered <- Y - Y.hat
  W.centered <- W - W.hat
  train.data <- create_train_matrices(X, outcome = Y.centered, treatment = W.centered)

  # Validate variables
  num.threads <- validate_num_threads(num.threads)
  linear.correction.variables <- validate_ll_vars(linear.correction.variables, ncol(X))
  ll.lambda <- validate_ll_path(lambda.path)

  # Subtract 1 to account for C++ indexing
  linear.correction.variables <- linear.correction.variables - 1

  args <- list(forest.object = forest.short,
               num.threads = num.threads,
               estimate.variance = FALSE,
               ll.lambda = ll.lambda,
               ll.weight.penalty = ll.weight.penalty,
               linear.correction.variables = linear.correction.variables)

  # Find sequence of predictions by lambda
  prediction.object <- do.call.rcpp(ll_causal_predict_oob, c(train.data, args))
  predictions <- prediction.object$predictions
  errors <- apply(predictions, MARGIN = 2, FUN = function(row) {
    # compute R-learner loss
    mean(((Y - Y.hat) - (W - W.hat) * row)**2)
  })

  return(list(
    lambdas = ll.lambda, errors = errors, oob.predictions = predictions,
    lambda.min = ll.lambda[which.min(errors)]
  ))
}
