#' Average LATE (deprecated)
#'
#' See the function `average_treatment_effect`
#'
#' @param forest The forest
#' @param ... Additional arguments (currently ignored).
#'
#' @return output
#'
#' @export
average_late <- function(forest, ...) {
  stop("This function has been deprecated after version 1.2.0. See the function `average_treatment_effect` instead.")
}

#' Average partial effect (deprecated)
#'
#' See the function `average_treatment_effect`
#'
#' @param forest The forest
#' @param ... Additional arguments (currently ignored).
#'
#' @return output
#'
#' @export
average_partial_effect <- function(forest, ...) {
  stop("This function has been deprecated after version 1.2.0. See the function `average_treatment_effect` instead.")
}

#' Regression forest tuning (deprecated)
#'
#' To tune a regression forest, see the function `regression_forest`
#'
#' @param X X
#' @param Y Y
#' @param ... Additional arguments (currently ignored).
#'
#' @return output
#'
#' @export
tune_regression_forest <- function(X, Y, ...) {
  stop("This function has been deprecated after version 1.2.0.")
}
