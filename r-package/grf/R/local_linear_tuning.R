#' Local linear forest tuning
#' 
#' Finds the optimal ridge penalty for local linear prediction.
#'
#' @param forest The forest used for prediction.
#' @return The optimal ridge penalty lambda.
#'
#' @examples \dontrun{
#' # Find the optimal tuning parameters.
#' n = 500; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' forest = regression_forest(X,Y)
#' tuned.lambda = tune_locally_linear_forest(forest)
#'
#' # Use this parameter to predict from a locally linear forest.
#' predictions = predict(forest, linear.correction.variables = 1:p, lambda = tuned.lambda)
#'
#' @export
tune_local_linear_forest <- function(forest, linear.correction.variables) {
  Y = forest[["Y.orig"]]
  X = forest[["X.orig"]]

  # validate correction variables and subtract 1 to account for C++ indexing
  linear.correction.variables = validate_vars(linear.correction.variables, ncol(X))
  linear.correction.variables = linear.correction.variables - 1

  vals = -10:5
  lambdas = exp(vals)

  errors = sapply(lambdas, function(lambda){
    oob.predictions = predict(forest, linear.correction.variables = linear.correction.variables, lambda = lambda)$predictions
    mean( (oob.predictions - Y)**2 )
  })
  tuned.lambda = lambdas[which.min(errors)]
  tuned.lambda
}
