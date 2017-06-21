#' Retrieve a single tree from a trained forest object.
#'
#' @param forest The trained forest.
#' @param index The index of the tree to retrieve.
#'
#' @return A GRF tree object.
#' @export
examine.tree = function(forest, index) {
	tree = examine_tree(forest, 500)
	class(tree) = "grf.tree"
	tree$columns = colnames(forest$original.data)
	tree
}
