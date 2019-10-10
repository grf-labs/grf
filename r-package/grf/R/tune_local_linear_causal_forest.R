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
#' \dontrun{
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
#' predictions <- predict(forest, linear.correction.variables = 1:p, lambda = tuned.lambda)
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

  data <- create_data_matrices(X, Y.centered, W.centered)

  # Validate variables
  num.threads <- validate_num_threads(num.threads)
  linear.correction.variables <- validate_ll_vars(linear.correction.variables, ncol(X))
  lambda.path <- validate_ll_path(lambda.path)

  # Subtract 1 to account for C++ indexing
  linear.correction.variables <- linear.correction.variables - 1

  # Find sequence of predictions by lambda
  prediction.object <- ll_causal_predict_oob(forest.short, data$train.matrix, data$sparse.train.matrix,
      data$outcome.index, data$treatment.index, lambda.path, ll.weight.penalty, linear.correction.variables, num.threads, FALSE)
  predictions <- prediction.object$predictions

  errors <- apply(predictions, MARGIN = 2, FUN = function(row) {
    # compute R-learner loss
    mean(((Y - Y.hat) - (W - W.hat) * row)**2)
  })

  return(list(
    lambdas = lambda.path, errors = errors, oob.predictions = predictions,
    lambda.min = lambda.path[which.min(errors)]
  ))
}
