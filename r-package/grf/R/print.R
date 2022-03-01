#' Print a GRF forest object.
#' @param x The tree to print.
#' @param decay.exponent A tuning parameter that controls the importance of split depth.
#' @param max.depth The maximum depth of splits to consider.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print grf
#' @export
print.grf <- function(x, decay.exponent = 2, max.depth = 4, ...) {
  var.importance <- variable_importance(x, decay.exponent, max.depth)
  var.importance <- c(round(var.importance, 3))
  names(var.importance) <- 1:length(var.importance)

  main.class <- class(x)[1]
  num.samples <- nrow(x$X.orig)

  cat("GRF forest object of type", main.class, "\n")
  cat("Number of trees:", x[["_num_trees"]], "\n")
  cat("Number of training samples:", num.samples, "\n")

  cat("Variable importance:", "\n")
  print(var.importance)
}

#' Print a GRF tree object.
#' @param x The tree to print.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print grf_tree
#' @export
print.grf_tree <- function(x, ...) {
  cat("GRF tree object", "\n")
  cat("Number of training samples:", x$num_samples, "\n")
  cat("Variable splits:", "\n")

  # Add the index of each node as an attribute for easy access.
  nodes <- lapply(1:length(x$nodes), function(i) {
    node <- x$nodes[[i]]
    node$index <- i
    return(node)
  })

  # Perform DFS to print the nodes (mimicking a stack with a list).
  frontier <- nodes[1]
  frontier[[1]]$depth <- 0
  while (length(frontier) > 0) {
    # Pop the first node off the stack.
    node <- frontier[[1]]
    frontier <- frontier[-1]

    output <- paste(rep("  ", node$depth), collapse = "")
    output <- paste(output, "(", node$index, ")", sep = "")

    if (node$is_leaf) {
      leaf_stats_text <- ""
      if (!is.null(node$leaf_stats)) {
        leaf_stats_text <- paste(paste(names(node$leaf_stats), unname(node$leaf_stats), sep = ": ", collapse = " "))
      }
      output <- paste(output, "* num_samples:", length(node$samples), "", leaf_stats_text)
    } else {
      split.var <- node$split_variable
      split.var.name <- x$columns[split.var]
      output <- paste(output, "split_variable:", split.var.name, " split_value:", signif(node$split_value))

      left_child <- nodes[node$left_child]
      left_child[[1]]$depth <- node$depth + 1

      right_child <- nodes[node$right_child]
      right_child[[1]]$depth <- node$depth + 1

      frontier <- c(left_child, right_child, frontier)
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
  cat("Number of forests: ", length(x$forests), "\n")
}


#' Print the Rank-Weighted Average Treatment Effect (RATE).
#' @param x The output of rank_average_treatment_effect.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print rank_average_treatment_effect
#' @export
print.rank_average_treatment_effect <- function(x, ...) {
  df <- data.frame(estimate = x[["estimate"]], std.err = x[["std.err"]], target = x[["target"]],
                   stringsAsFactors = FALSE)
  print(df, row.names = FALSE)
  invisible(x)
}


#' Print tuning output.
#' Displays average error for q-quantiles of tuned parameters.
#' @param x The tuning output to print.
#' @param tuning.quantiles vector of quantiles to display average error over.
#'  Default: seq(0, 1, 0.2) (quintiles)
#' @param ... Additional arguments (currently ignored).
#'
#' @method print tuning_output
#' @export
print.tuning_output <- function(x, tuning.quantiles = seq(0, 1, 0.2), ...) {
  if (x$status == "failure") {
    cat("Tuning status: failure.\n")
    cat("This indicates tuning was attempted but failed due to an error, and we fell back to default parameters: \n\n")
    params <- x$params
    cat(paste0(names(params), ": ", params, "\n"))
  } else if (x$status == "default") {
    cat("Tuning status: default.\n")
    cat("This indicates tuning was attempted. ")
    cat("However, we could not find parameters that were expected to perform better than default: \n\n")
    params <- x$params
    cat(paste0(names(params), ": ", params, "\n"))
  } else if (x$status == "tuned") {
    cat("Tuning status: tuned.\n")
    cat("This indicates tuning found parameters that are expected to perform better than default. \n\n")
    grid <- x$grid
    out <- lapply(colnames(grid)[-1], function(name) {
      q <- stats::quantile(grid[, name], probs = tuning.quantiles)
      # Cannot form for example quintiles for mtry if the number of variables is
      # less than 5, so here we just truncate the groups.
      if (length(unique(q) < length(q))) {
        q <- unique(q)
      }
      rank <- cut(grid[, name], q, include.lowest = TRUE)
      # If the cut is only a single level, the variable is binary and we can
      # aggregate by that directly.
      if (length(levels(rank)) == 1) {
        rank = grid[, name]
      }
      out <- stats::aggregate(grid[, "error"], by = list(rank), FUN = mean)
      colnames(out) <- c(name, "error")
      out
    })

    cat(paste0("Predicted debiased error: ", x$error, "\n\n"))

    cat("Tuned parameters: \n")
    cat(paste0(names(x$params), ": ", x$params, "\n"))
    cat("\n")

    cat("Average error by ", length(tuning.quantiles) - 1, "-quantile:\n", sep = "")
    for (i in out) {
      cat("\n")
      print(i, row.names = FALSE)
    }
  } else {
    stop(paste0("Error while reading tuning output. ",
                "Parameter 'status' must be one of 'failure', 'default', or 'tuned'"))
  }
}
