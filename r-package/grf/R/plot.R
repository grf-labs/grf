#' Writes each node information
#' If it is a leaf node: show it in different color, show number of samples, show leaf id
#' If it is a non-leaf node: show its splitting variable and splitting value
#' @param tree the tree to convert
#' @param index the index of the current node
#' @keyword internal
create_dot_body <- function(tree, index = 1) {
  node <- tree$nodes[[index]]

  # Leaf case: print label only
  if (node$is_leaf) {
    num_samples <- length(node$samples)
    leaf_stats_text <- ""
    if(!is.null(node$leaf_stats)){
      leaf_stats_text <- paste("\n", paste(node$leaf_stats$label, node$leaf_stats$value, sep = " = ", collapse = "\n"))
    }
    line_label <- paste(index - 1, ' [shape=box,style=filled,color=".7 .3 1.0" , label="size = ',
                        num_samples, leaf_stats_text, '"];')
    return(line_label)
  }

  # Non-leaf case: print label, child edges
  if (!is.null(node$left_child)) {
    edge <- paste(index - 1, "->", node$left_child - 1)
    if (index == 1) {
      edge_info_left <- paste(edge, '[labeldistance=2.5, labelangle=45, headlabel="True"];')
    }
    else {
      edge_info_left <- paste(edge, " ;")
    }
  }
  else {
    edge_info_right <- NULL
  }

  if (!is.null(node$right_child)) {
    edge <- paste(index - 1, "->", node$right_child - 1)
    if (index == 1) {
      edge_info_right <- paste(edge, '[labeldistance=2.5, labelangle=-45, headlabel="False"]')
    } else {
      edge_info_right <- paste(edge, " ;")
    }
  } else {
    edge_info_right <- NULL
  }

  variable_name <- tree$columns[node$split_variable]
  node_info <- paste(index - 1, '[label="', variable_name, "<=", round(node$split_value, 2), '"] ;')

  this_lines <- paste(node_info,
    edge_info_left,
    edge_info_right,
    sep = "\n"
  )

  left_child_lines <- ifelse(!is.null(node$left_child),
    create_dot_body(tree, index = node$left_child),
    NULL
  )

  right_child_lines <- ifelse(!is.null(node$right_child),
    create_dot_body(tree, index = node$right_child),
    NULL
  )

  lines <- paste(this_lines, left_child_lines, right_child_lines, sep = "\n")

  return(lines)
}

#' Export a tree in DOT format.
#' This function generates a GraphViz representation of the tree,
#' which is then written into `dot_string`.
#' @param tree the tree to convert
#' @keyword internal
export_graphviz <- function(tree) {
  header <- "digraph nodes { \n node [shape=box] ;"
  footer <- "}"
  body <- create_dot_body(tree)

  dot_string <- paste(header, body, footer, sep = "\n")

  return(dot_string)
}

#' Plot a GRF tree object.
#' @param x The tree to plot
#' @param ... Additional arguments (currently ignored).
#'
#' @method plot grf_tree
#' @export
plot.grf_tree <- function(x, ...) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package \"DiagrammeR\" must be installed to plot trees.")
  }

  dot_file <- export_graphviz(x)
  DiagrammeR::grViz(dot_file)
}

leaf_stats <- function(forest, samples) UseMethod("leaf_stats")


#' A default leaf_stats for forests classes without a leaf_stats method
#' that always returns NULL.
#' @param forest Any forest
#' @param samples The samples to include in the calculations.
#' @param ... Additional arguments (currently ignored).
#'
#' @return NULL
#'
#' @method leaf_stats default
leaf_stats.default <- function(forest, samples, ...){
  return(NULL)
}

#' Calculate summary stats given a set of samples for quantile forests.
#' @param forest The GRF forest
#' @param samples The samples to include in the calculations.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A label, value dataframe containing summary stats
#'
#' @method leaf_stats quantile_forest
leaf_stats.quantile_forest <- function(forest, samples, ...){
  funcs <- c(
    function(forest, samples){
      label <- "average_X"
      res <- round(mean(forest$X.orig[samples]), 2)
      return(data.frame(label = label, value = res))
    },
    function(forest, samples){
      label <- "average_Y"
      res <- round(mean(forest$Y.orig[samples]), 2)
      return(data.frame(label = label, value = res))
    }
  )
  return(calc_leaf_stats(forest, samples, funcs))
}

#' Calculate summary stats given a set of samples for causal forests.
#' @param forest The GRF forest
#' @param samples The samples to include in the calculations.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A label, value dataframe containing summary stats
#'
#' @method leaf_stats causal_forest
leaf_stats.causal_forest <- function(forest, samples, ...){
  funcs <- c(
    function(forest, samples){
      label <- "average_W"
      res <- round(mean(forest$W.orig[samples]), 2)
      return(data.frame(label = label, value = res))
    },
    function(forest, samples){
      label <- "average_W_hat"
      res <- round(mean(forest$W.hat[samples]), 2)
      return(data.frame(label = label, value = res))
    },
    function(forest, samples){
      label <- "average_Y"
      res <- round(mean(forest$Y.orig[samples]), 2)
      return(data.frame(label = label, value = res))
    },
    function(forest, samples){
      label <- "average_Y_hat"
      res <- round(mean(forest$Y.hat[samples]), 2)
      return(data.frame(label = label, value = res))
    }
  )
  return(calc_leaf_stats(forest, samples, funcs))
}

#' Helper function to calculate an arbitrary set of summary stats.
#' @param forest The GRF forest
#' @param samples The samples to include in the calculations
#' @param funcs A vector of functions that return label,value dataframes
#' @param ... Additional arguments (currently ignored).
#'
#' @return A label, value dataframe containing summary stats
calc_leaf_stats <- function(forest, samples, funcs){
  res <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("label", "value"))
  for(func in funcs){
    res <- rbind(res, func(forest, samples))
  }
  return(res)
}
