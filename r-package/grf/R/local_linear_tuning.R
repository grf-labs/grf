#' Local linear forest tuning
#' 
#' Finds the optimal ridge penalty for local linear prediction.
#'
#' @param forest The forest used for prediction.
#' @param linear.correction.variables Variables to use for local linear prediction. If left null,
#'          all variables are used.
#' @param ll.weight.penalty Option to standardize ridge penalty by covariance (TRUE),
#'                            or penalize all covariates equally (FALSE). Defaults to FALSE.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param lambda.path Optional list of lambdas to use for cross-validation.
#' @return A list of lambdas tried, corresponding errors, and optimal ridge penalty lambda.
#'
#' @examples \dontrun{
#' # Find the optimal tuning parameters.
#' n = 500; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' forest = regression_forest(X,Y)
#' tuned.lambda = tune_ll_regression_forest(forest)
#'
#' # Use this parameter to predict from a local linear forest.
#' predictions = predict(forest, linear.correction.variables = 1:p, lambda = tuned.lambda)
#' }
#'
#' @export
tune_ll_regression_forest <- function(forest,
                                     linear.correction.variables = NULL,
                                     ll.weight.penalty = FALSE,
                                     num.threads = NULL,
                                     lambda.path = NULL) {
  forest.short = forest[-which(names(forest) == "X.orig")]

  X = forest[["X.orig"]]
  Y = forest[["Y.orig"]]
  data = create_data_matrices(X, Y)
  outcome.index = ncol(X) + 1

  # Validate variables
  num.threads = validate_num_threads(num.threads)
  linear.correction.variables = validate_ll_vars(linear.correction.variables, ncol(X))
  lambda.path = validate_ll_path(lambda.path)

  # Subtract 1 to account for C++ indexing
  linear.correction.variables = linear.correction.variables - 1

  # Enforce no variance estimates in tuning
  estimate.variance = FALSE
  prediction.object = ll_regression_predict_oob(forest.short, data$default, data$sparse, outcome.index,
      lambda.path, ll.weight.penalty, linear.correction.variables, num.threads, estimate.variance)

  prediction.object = prediction.object$predictions
  errors = apply(prediction.object, MARGIN = 2, FUN = function(row){
    mean( (row - Y)**2 )
  })

  return(list(lambdas = lambda.path, errors = errors, oob.predictions = prediction.object,
           lambda.min = lambda.path[which.min(errors)]))
}
