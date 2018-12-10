#' Merges a list of forests that were grown using the same data into one large forest.
#'
#' @param forest_list A `list` of forests to be concatenated. 
#'                        All forests must be of the same type, and the type must be a subclass of `grf`.
#'                        
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed.
#'        Note that even if OOB predictions have already been precomputed for the forests in 'forest_list',
#'        those predictions are not used. Instead, a new set of oob predictions is computed anew using the 
#'        larger forest.
#' 
#' @return A single forest containing all the trees in each forest in the input list.
#'
#' @examples \dontrun{
#' # Train standard regression forests
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.forest1 = regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 100)
#' r.forest2 = regression_forest(X, Y, compute.oob.predictions = FALSE, num.trees = 100)
#'
#' # Join the forests together. The resulting forest will contain 200 trees.
#' big_rf = merge_forests(list(r.forest1, r.forest2))
#' }
#' 
#' @export
merge_forests <- function(forest_list, compute.oob.predictions=TRUE) {
  
  validate_forest_list(forest_list)
  
  first_forest <- forest_list[[1]]
  data <- create_data_matrices(first_forest$X.orig, first_forest$Y.orig)
  
  big_forest <- merge(forest_list, data$default, data$sparse)
  
  big_forest[["ci.group.size"]] <- first_forest$ci.group.size
  big_forest[["X.orig"]] <- first_forest$X.orig
  big_forest[["Y.orig"]] <- first_forest$Y.orig
  big_forest[["clusters"]] <- first_forest$clusters
  big_forest[["min.node.size"]] <- first_forest$min.node.size
  
  class(big_forest) <- class(first_forest)
  
  if (compute.oob.predictions) {
<<<<<<< HEAD:r-package/grf/R/merge_forests.R
    oob.pred <- predict(big_forest)
    big_forest[["predictions"]] <- oob.pred$predictions
    # Must include checks here because some big_forest types may have this 
    # method not yet implemented
    if (!is.null(big_forest["debiased.error"])) {
      big_forest[["debiased.error"]] <- oob.pred$debiased.error
    }
    if (!is.null(big_forest["excess.error"])) {
      big_forest[["excess.error"]] <- oob.pred$excess.error
=======
    oob.pred <- predict(forest)
    forest[["predictions"]] <- oob.pred$predictions
    # Must include checks here because some forest types may have this 
    # method not yet implemented
    if (!is.null(forest["debiased.error"])) {
      forest[["debiased.error"]] <- oob.pred$debiased.error
    }
    if (!is.null(forest["excess.error"])) {
      forest[["excess.error"]] <- oob.pred$excess.error
>>>>>>> dcde786b568f06da3194c5122242b0ea09aa9de9:r-package/grf/R/merge_forests.R
    }
  }
  
  big_forest
}


validate_forest_list <- function(forest_list) {
  
  if (length(forest_list) == 0) {
    stop("Length of argument 'forest_list' must be positive.")
  }
  
  first_forest <- forest_list[[1]]
  if (!is(first_forest, "grf")) {
    stop("Argument 'forest_list' must be a list of grf objects. 
           Be sure to use 'list(forest1, forest2), not 'c(forest1, forest2)'.")
  }
  
  classes <- unique(sapply(forest_list, class)[1,])
  if (length(classes) > 1) {
    stop(paste("All forests in 'forest_list' must be of the same type, but we found:", 
               paste(classes, collapse=", ")))
  }
  
}