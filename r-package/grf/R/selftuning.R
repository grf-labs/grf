#' Selftuning
#' 
#' selftuning parameters for random forest
#' parameters including
#' mrty
#' min_node_size
#' sample_fraction
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' 
#' @return A trained regression forest object.
#'
#' @examples
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
#'
#' @export
selftuning <- function(X, Y, mtry_l = 0.01, mtry_r = 1, min_node_size_l = 1, min_node_size_r = 100, sample_faction_l = 0.01, sample_faction_r = 0.5) {
  
  if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
  
  input.data <- as.matrix(cbind(X, Y))
  variable.names <- c(colnames(X), "outcome")
  outcome.index <- ncol(input.data)
  
  bestPara <- selftuning_train(input.data, outcome.index, variable.names, mtry_l, mtry_r, min_node_size_l, min_node_size_r, sample_faction_l, sample_faction_r)
}