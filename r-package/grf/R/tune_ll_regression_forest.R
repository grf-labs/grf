#' Local linear forest tuning
#'
#' Finds the optimal ridge penalty for local linear prediction.
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
#' @keywords internal
tune_ll_regression_forest <- function(forest,
                                      linear.correction.variables = NULL,
                                      ll.weight.penalty = FALSE,
                                      num.threads = NULL,
                                      lambda.path = NULL) {
  forest.short <- forest[-which(names(forest) == "X.orig")]
  X <- forest[["X.orig"]]
  Y <- forest[["Y.orig"]]
  train.data <- create_train_matrices(X, outcome = Y)

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

  prediction.object <- do.call.rcpp(ll_regression_predict_oob, c(train.data, args))
  predictions <- prediction.object$predictions
  errors <- apply(predictions, MARGIN = 2, FUN = function(row) {
    mean((row - Y)**2)
  })

  return(list(
    lambdas = ll.lambda, errors = errors, oob.predictions = predictions,
    lambda.min = ll.lambda[which.min(errors)]
  ))
}
