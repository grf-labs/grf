#' Retrieve a single tree from a trained forest object.
#'
#' @param forest The trained forest.
#' @param index The index of the tree to retrieve.
#'
#' @return A GRF tree object.
#' @export
get_tree = function(forest, index) {
	if (index < 1 || index > forest$num.trees) {
		stop(paste("The provided index,", index, "is not valid."))
	}

	tree = deserialize_tree(forest, 500)
	class(tree) = "grf_tree"
	tree$columns = colnames(forest$original.data)
	tree
}

#' Get summaries of which features the forest split on
#'
#' @param forest The trained forest.
#' @param max.depth Maximum depth of splits to consider.
#'
#' @export
split_frequencies = function(forest, max.depth=4) {
  raw = compute_split_frequencies(forest, max.depth)
  raw[,forest$feature.indices]
}
