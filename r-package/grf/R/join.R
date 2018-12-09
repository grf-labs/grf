#' Concatenates a list of identical forests that were grown using the same data. 
#' Useful for training forests in a distributed manner.
#'
#' @param regression_list A `list` of forests to be concatenated. 
#'                        All forests must be of the same type, and the type must be a subclass of `grf`.
#'                        
#' 
#' @return A larger forest that is the concatenation of all previous forests.
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
#' big_rf = join_forests(list(r.forest1, r.forest2))
#' }
#' 
#' @export
join_forests <- function(forest_list) {
  
  validate_forest_list(forest_list)
  
  first_forest <- forest_list[[1]]
  data <- create_data_matrices(first_forest$X.orig, first_forest$Y.orig)
  
  big_forest <- cpp_join_forests(forest_list, data$default, data$sparse)
  
  big_forest[["ci.group.size"]] <- first_forest$ci.group.size
  big_forest[["X.orig"]] <- first_forest$X.orig
  big_forest[["Y.orig"]] <- first_forest$Y.orig
  big_forest[["clusters"]] <- first_forest$clusters
  big_forest[["min.node.size"]] <- first_forest$min.node.size
  
  class(big_forest) <- class(first_forest)
  
  big_forest
}


validate_forest_list <- function(forest_list) {
  
  if (length(forest_list) == 0) {
    stop("Length of argument 'forest_list' must be positive.")
  }
  
  first_forest <- forest_list[[1]]
  if (!is(rf1, "grf")) {
    stop("Argument 'forest_list' must be a list of grf objects. 
           Be sure to use 'list(forest1, forest2), not 'c(forest1, forest2)'.")
  }
  
  classes <- unique(sapply(forest_list, class)[1,])
  if (length(classes) > 1) {
    stop(paste("All forests in 'forest_list' must be of the same type, but we found:", 
               paste(classes, collapse=", ")))
  }
  
}