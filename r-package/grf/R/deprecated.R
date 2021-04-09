#' Average LATE (removed)
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
  stop("This function has been removed after version 1.2.0. See the function `average_treatment_effect` instead.")
}

#' Average partial effect (removed)
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
  stop("This function has been removed after version 1.2.0. See the function `average_treatment_effect` instead.")
}

#' Regression forest tuning (removed)
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
  stop("This function has been removed after version 1.2.0. To tune a regression_forest see the `tune.parameters` argument in that forest.")
}

#' Causal forest tuning (removed)
#'
#' To tune a causal forest, see the function `causal_forest`
#'
#' @param X X
#' @param Y Y
#' @param W W
#' @param Y.hat Y.hat
#' @param W.hat W.hat
#' @param ... Additional arguments (currently ignored).
#'
#' @return output
#'
#' @export
tune_causal_forest <- function(X, Y, W, Y.hat, W.hat, ...) {
  stop("This function has been removed after version 1.2.0. To tune a causal_forest see the `tune.parameters` argument in that forest.")
}

#' Instrumental forest tuning (removed)
#'
#' To tune a instrumental forest, see the function `instrumental_forest`
#'
#' @param X X
#' @param Y Y
#' @param W W
#' @param Z Z
#' @param Y.hat Y.hat
#' @param W.hat W.hat
#' @param Z.hat Z.hat
#' @param ... Additional arguments (currently ignored).
#'
#' @return output
#'
#' @export
tune_instrumental_forest <- function(X, Y, W, Z, Y.hat, W.hat, Z.hat, ...) {
  stop("This function has been removed after version 1.2.0. To tune an instrumental_forest see the `tune.parameters` argument in that forest.")
}

#' Custom forest (removed)
#'
#' To build a custom forest, see an existing simpler forest, like regression_forest,
#' for a development template.
#'
#' @param X X
#' @param Y Y
#' @param ... Additional arguments (currently ignored).
#'
#' @export
custom_forest <- function(X, Y, ...) {
  stop("This function has been removed after version 1.2.0. For a development template, see existing forests, such as regression_forest, for inspiration.")
}

#' Retrieve forest weights (renamed to get_forest_weights)
#'
#' @param forest The trained forest.
#' @param ... Additional arguments (currently ignored).
#'
#' @export
get_sample_weights <- function(forest, ...) {
  stop("This function has been renamed after version 1.2.0. See the function `get_forest_weights`.")
}
