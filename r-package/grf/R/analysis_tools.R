#' Retrieve a single tree from a trained forest object.
#'
#' @param forest The trained forest.
#' @param index The index of the tree to retrieve.
#'
#' @return A GRF tree object.
#'
#' @examples \dontrun{
#' # Train a quantile forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
#'
#' # Examine a particular tree.
#' q.tree = get_tree(q.forest, 3)
#' q.tree$nodes
#' }
#'
#' @export
get_tree = function(forest, index) {
	if (index < 1 || index > forest$num.trees) {
		stop(paste("The provided index,", index, "is not valid."))
	}

	tree = deserialize_tree(forest, index)
	class(tree) = "grf_tree"

	columns = colnames(forest$X.orig)
	indices = 1:ncol(forest$X.orig)
	tree$columns  = sapply(indices, function(i) {
		if (!is.null(columns) & length(columns[i]) > 0) columns[i]
		else paste("X", i, sep=".")
	})

	tree
}

#' Calculate which features the forest split on at each depth.
#'
#' @param forest The trained forest.
#' @param max.depth Maximum depth of splits to consider.
#'
#' @return A matrix of split depth by feature index, where each value
#' is the number of times the feature was split on at that depth.
#'
#' @examples \dontrun{
#' # Train a quantile forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
#'
#' # Calculate the split frequencies for this forest.
#' split_frequencies(q.forest)
#' }
#'
#' @export
split_frequencies = function(forest, max.depth=4) {
  raw = compute_split_frequencies(forest, max.depth)
  feature.indices = 1:ncol(forest$X.orig)
  raw[,feature.indices, drop = FALSE]
}

#' Calculate a simple measure of 'importance' for each feature.
#'
#' @param forest The trained forest.
#' @param decay.exponent A tuning parameter that controls the importance of split depth.
#' @param max.depth Maximum depth of splits to consider.
#'
#' @return A list specifying an 'importance value' for each feature.
#'
#' @examples \dontrun{
#' # Train a quantile forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
#'
#' # Calculate the 'importance' of each feature.
#' variable_importance(q.forest)
#' }
#'
#' @export
variable_importance = function(forest, decay.exponent=2, max.depth=4) {
  split.freq <- split_frequencies(forest, max.depth)
  split.freq <- split.freq / pmax(1L, rowSums(split.freq))
  weight <- seq_len(nrow(split.freq)) ^ -decay.exponent
  t(split.freq) %*% weight / sum(weight)
}
