#' Retrieve a single tree from a trained forest object.
#'
#' @param forest The trained forest.
#' @param index The index of the tree to retrieve.
#'
#' @return A GRF tree object.
#'
#' @examples
#' # Train a quantile forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
#'
#' # Examine a particular tree.
#' q.tree = get_tree(q.forest, 3)
#' q.tree$nodes
#'
#' @export
get_tree = function(forest, index) {
	if (index < 1 || index > forest$num.trees) {
		stop(paste("The provided index,", index, "is not valid."))
	}

	tree = deserialize_tree(forest, index)
	class(tree) = "grf_tree"
	tree$columns = colnames(forest$original.data)
	tree
}

#' Get summaries of which features the forest split on
#'
#' @param forest The trained forest.
#' @param max.depth Maximum depth of splits to consider.
#'
#' @return A matrix of split depth by feature index, where each value
#' is the number of times the feature was split on at that depth.
#'
#' @examples
#' # Train a quantile forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
#'
#' # Calculate the split frequencies for this forest.
#' split_frequencies(q.forest)
#'
#' @export
split_frequencies = function(forest, max.depth=4) {
  raw = compute_split_frequencies(forest, max.depth)
  raw[,forest$feature.indices]
}
