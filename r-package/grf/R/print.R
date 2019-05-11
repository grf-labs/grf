#' Print a GRF forest object.
#' @param x The tree to print.
#' @param decay.exponent A tuning parameter that controls the importance of split depth.
#' @param max.depth The maximum depth of splits to consider.
#' @param tuning.quantiles Number of quantiles to display average error over for tuned parameters.
#'                         Default: 0 (no output).
#' @param ... Additional arguments (currently ignored).
#'
#' @method print grf
#' @export
print.grf <- function(x, decay.exponent=2, max.depth=4, tuning.quantiles=0, ...) {
    var.importance = variable_importance(x, decay.exponent, max.depth)
    var.importance = c(round(var.importance, 3))
    names(var.importance) = 1:length(var.importance)

    main.class = class(x)[1]
    num.samples= nrow(x$X.orig)

    cat("GRF forest object of type", main.class, "\n")
    cat("Number of trees: ", x[["_num_trees"]], "\n")
    cat("Number of training samples:", num.samples, "\n")

    cat("Variable importance:", "\n")
    print(var.importance)

    if (exists('tuning.output', x)) {
      if (tuning.quantiles > 0) {
        p = ncol(x$X.orig)
        print_tuning_params(x$tuning.output, tuning.quantiles, p)
      }
    }
}

#' Print a GRF tree object.
#' @param x The tree to print.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print grf_tree
#' @export
print.grf_tree <- function(x, ...) {
    cat("GRF tree object", "\n")
    cat("Number of training samples: ", x$num_samples, "\n")
    cat("Variable splits:", "\n")

    # Add the index of each node as an attribute for easy access.
    nodes = lapply(1:length(x$nodes), function(i) {
        node = x$nodes[[i]]
        node$index = i
        return(node)
    })

    # Perform DFS to print the nodes (mimicking a stack with a list).
    frontier = nodes[1]
    frontier[[1]]$depth = 0
    while (length(frontier) > 0) {
        # Pop the first node off the stack.
        node = frontier[[1]]
        frontier = frontier[-1]

        output = paste(rep("  ", node$depth), collapse="")
        output = paste(output, "(", node$index, ")", sep="")

        if (node$is_leaf) {
            output = paste(output, "* num_samples:", length(node$samples))
        } else {
            split.var = node$split_variable
            split.var.name = x$columns[split.var]
            output = paste(output, "split_variable:", split.var.name, " split_value:", signif(node$split_value))

            left_child = nodes[node$left_child]
            left_child[[1]]$depth = node$depth + 1

            right_child = nodes[node$right_child]
            right_child[[1]]$depth = node$depth + 1

            frontier = c(left_child, right_child, frontier)
        }
        cat(output, "\n")
    }
}


#' Print a boosted regression forest
#' @param x The boosted forest to print.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print boosted_regression_forest
#' @export
print.boosted_regression_forest <- function(x, ...) {
    cat("Boosted GRF object", "\n")
    cat("Number of forests: ",length(x$forests), "\n")
}


print_tuning_params <- function(tuning.output, nq, p) {
  # Print average error for `nq`-quantiles of tuned parameters
  # Extra branches for `mtry` and `min.node.size`
  grid = tuning.output$grid

  out = lapply(colnames(grid)[-1], function(name) {
    m = nq
    if (name == "mtry")
      m = min(p - 1, nq)
    probs = 0:m/m
    q = quantile(grid[, name], probs = probs)
    if (length(unique(q) < length(q)))
      q = unique(q)
    rank = cut(grid[, name], q, include.lowest=TRUE)
    out = aggregate(grid[, "error"], by=list(rank), FUN=mean)
    colnames(out) = c(name, "error")
    out
  })

  err = tuning.output$error
  params = tuning.output$params[colnames(grid)[-1]]
  opt = formatC(c(err, params))

  cat("Optimal tuning parameters: \n")
  cat(paste0(names(opt), ": ", opt, "\n"))

  cat("Average error by ", nq, "-quantile:\n", sep="")
  for (i in out) {
    print(i, row.names = FALSE)
  }
}
