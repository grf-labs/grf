#' Merges a list of forests that were grown using the same data into one large forest.
#'
#' @param forest_list A `list` of forests to be concatenated.
#'                        All forests must be of the same type, and the type must be a subclass of `grf`.
#'                        In addition, all forests must have the same 'ci.group.size'.
#'                        Other tuning parameters (e.g. alpha, mtry, min.node.size, imbalance.penalty) are
#'                        allowed to differ across forests.
#'
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed.
#'        Note that even if OOB predictions have already been precomputed for the forests in 'forest_list',
#'        those predictions are not used. Instead, a new set of oob predictions is computed anew using the
#'        larger forest. Default is TRUE.
#'
#' @return A single forest containing all the trees in each forest in the input list.
#'
#' @examples
#' \dontrun{
#' # Train standard regression forests
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' r.forest1 <- regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 100)
#' r.forest2 <- regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 100)
#'
#' # Join the forests together. The resulting forest will contain 200 trees.
#' big_rf <- merge_forests(list(r.forest1, r.forest2))
#' }
#'
#' @export
merge_forests <- function(forest_list, compute.oob.predictions = TRUE) {
  validate_forest_list(forest_list)
  first_forest <- forest_list[[1]]

  big_forest <- merge(forest_list)

  # Make sure the new forest contains the necessary saved components like 'X.orig'.
  class(big_forest) <- class(first_forest)
  for (name in names(first_forest)) {
    if (!startsWith(name, "_")
    && name != "predictions"
    && name != "debiased.error"
    && name != "excess.error") {
      big_forest[[name]] <- first_forest[[name]]
    }
  }

  if (compute.oob.predictions) {
    oob.pred <- predict(big_forest)
    big_forest[["predictions"]] <- oob.pred$predictions
    big_forest[["debiased.error"]] <- oob.pred$debiased.error
    big_forest[["excess.error"]] <- oob.pred$excess.error
  }

  big_forest
}

#' @importFrom methods is
validate_forest_list <- function(forest_list) {
  if (length(forest_list) == 0) {
    stop("Length of argument 'forest_list' must be positive.")
  }

  first_forest <- forest_list[[1]]
  if (!is(first_forest, "grf")) {
    stop("Argument 'forest_list' must be a list of grf objects.
           Be sure to use 'list(forest1, forest2), not 'c(forest1, forest2)'.")
  }

  classes <- unique(sapply(forest_list, class)[1, ])
  if (length(classes) > 1) {
    stop(paste(
      "All forests in 'forest_list' must be of the same type, but we found:",
      paste(classes, collapse = ", ")
    ))
  }

  n.cols <- unique(lapply(forest_list, function(x) {ncol(x$X.orig)}))
  n.obs <- unique(lapply(forest_list, function(x) {nrow(x$X.orig)}))
  if (length(n.cols) != 1 || length(n.obs) != 1) {
    stop("All forests in 'forest_list' must be trained on the same data.")
  }
}
