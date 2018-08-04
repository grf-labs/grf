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

  vals = -10:5
  lambdas = exp(vals)

  preds.matrix = predict(forest, linear.correction.variables = linear.correction.variables, lambda = lambdas)$predictions
  errors = sapply(1:length(lambdas), function(i){
    preds.lambda = preds.matrix[,i]
    mean( (preds.lambda - Y)**2 )
  })

  tuned.lambda = lambdas[which.min(errors)]
  c(tuned.lambda)
}
