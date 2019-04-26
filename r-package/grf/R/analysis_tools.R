#' Retrieve a single tree from a trained forest object.
#'
#' @param forest The trained forest.
#' @param index The index of the tree to retrieve.
#'
#' @return A GRF tree object containing the below attributes.
#'     drawn_samples: a list of examples that were used in training the tree. This includes
#'     examples that were used in choosing splits, as well as the examples that populate the leaf
#'     nodes. Put another way, if honesty is enabled, this list includes both subsamples from the
#'     split (J1 and J2 in the notation of the paper).
#'     num_samples: the number of examples used in training the tree.
#'     nodes: a list of objects representing the nodes in the tree, starting with the root node. Each
#'     node will contain an 'is_leaf' attribute, which indicates whether it is an interior or leaf node.
#'     Interior nodes contain the attributes 'left_child' and 'right_child', which give the indices of
#'     their children in the list, as well as 'split_variable', and 'split_value', which describe the
#'     split that was chosen. Leaf nodes only have the attribute 'samples', which is a list of the
#'     training examples that the leaf contains. Note that if honesty is enabled, this list will only
#'     contain examples from the second subsample that was used to 'repopulate' the tree (J2 in the
#'     notation of the paper).
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
	if (index < 1 || index > forest[["_num_trees"]]) {
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

#' Given a trained forest and test data, compute the training sample weights for each test point.
#'
#' During normal prediction, these weights are computed as an intermediate step towards producing estimates.
#' This function allows for examining the weights directly, so they could be potentially be used as the
#' input to a different analysis.
#'
#' @param forest The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).#' @param max.depth Maximum depth of splits to consider.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @return A sparse matrix where each row represents a test sample, and each column is a sample in the
#'         training data. The value at (i, j) gives the weight of training sample j for test sample i.
#'
#' @examples \dontrun{
#'  p = 10
#'  n = 100
#'  X = matrix(2 * runif(n * p) - 1, n, p)
#'  Y = (X[,1] > 0) + 2 * rnorm(n)
#'  rrf = regression_forest(X, Y, mtry=p)
#'  sample.weights.oob = get_sample_weights(rrf)
#'
#'  n.test = 15
#'  X.test = matrix(2 * runif(n.test * p) - 1, n.test, p)
#'  sample.weights = get_sample_weights(rrf, X.test)
#' }
#'
#' @export
get_sample_weights = function(forest, newdata = NULL, num.threads=NULL) {
  num.threads <- validate_num_threads(num.threads)

  forest.short <- forest[-which(names(forest) == "X.orig")]
  train.data <- create_data_matrices(forest[["X.orig"]])
  
  if (!is.null(newdata)) {
    data <- create_data_matrices(newdata)
    compute_weights(forest.short, train.data$default, train.data$sparse,
        data$default, data$sparse, num.threads)
  } else {
    compute_weights_oob(forest.short, train.data$default, train.data$sparse, num.threads)
  }
}
