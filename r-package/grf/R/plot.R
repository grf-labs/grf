#' Writes each node information
#' If it is a leaf node: show it in different color, show number of samples, show leaf id
#' If it is a non-leaf node: show its splitting variable and splitting value
#' If trained with missing values, the edge arrow is filled according to which direction the NAs are sent.
#' @param tree the tree to convert
#' @param index the index of the current node
#' @param include.na.path A boolean toggling whether to include the path of missing values or not.
#' @keywords internal
create_dot_body <- function(tree, index = 1, include.na.path) {
  node <- tree$nodes[[index]]

  # Leaf case: print label only
  if (node$is_leaf) {
    num_samples <- length(node$samples)
    leaf_stats_text <- ""
    if(!is.null(node$leaf_stats)){
      leaf_stats_text <- paste("\n", paste(names(node$leaf_stats), unname(node$leaf_stats), sep = " = ", collapse = "\n"))
    }
    line_label <- paste(index - 1, ' [shape=box,style=filled,color=".7 .3 1.0" , label="size = ',
                        num_samples, leaf_stats_text, '"];')
    return(line_label)
  }

  # Non-leaf case: print label, child edges
  if (!is.null(node$left_child)) {
    if (include.na.path) {
      na.left <- node$send_missing_left
      arrowhead <- ifelse(na.left, 'arrowhead=empty];', 'arrowhead=normal];')
    } else {
      arrowhead <- 'arrowhead=normal];'
    }
    edge <- paste(index - 1, "->", node$left_child - 1)
    if (index == 1) {
      edge_info_left <- paste(edge, '[labeldistance=2.5, labelangle=45, headlabel="True",')
    }
    else {
      edge_info_left <- paste(edge, '[')
    }
    edge_info_left <- paste(edge_info_left, arrowhead)
  }
  else {
    edge_info_right <- NULL
  }

  if (!is.null(node$right_child)) {
    if (include.na.path) {
      na.left <- node$send_missing_left
      arrowhead <- ifelse(na.left, 'arrowhead=normal];', 'arrowhead=empty];')
    } else {
      arrowhead <- 'arrowhead=normal];'
    }
    edge <- paste(index - 1, "->", node$right_child - 1)
    if (index == 1) {
      edge_info_right <- paste(edge, '[labeldistance=2.5, labelangle=-45, headlabel="False",')
    } else {
      edge_info_right <- paste(edge, '[')
    }
    edge_info_right <- paste(edge_info_right, arrowhead)
  } else {
    edge_info_right <- NULL
  }

  variable_name <- tree$columns[node$split_variable]
  split_value <- node$split_value
  node_info <- paste(index - 1, '[label="', variable_name, ifelse(is.na(split_value),"=", "<="), round(split_value, 2), '"] ;')

  this_lines <- paste(node_info,
    edge_info_left,
    edge_info_right,
    sep = "\n"
  )

  left_child_lines <- ifelse(!is.null(node$left_child),
    create_dot_body(tree, index = node$left_child, include.na.path = include.na.path),
    NULL
  )

  right_child_lines <- ifelse(!is.null(node$right_child),
    create_dot_body(tree, index = node$right_child, include.na.path = include.na.path),
    NULL
  )

  lines <- paste(this_lines, left_child_lines, right_child_lines, sep = "\n")

  return(lines)
}

#' Export a tree in DOT format.
#' This function generates a GraphViz representation of the tree,
#' which is then written into `dot_string`.
#' @param tree the tree to convert
#' @param include.na.path A boolean toggling whether to include the path of missing values or not.
#' @keywords internal
export_graphviz <- function(tree, include.na.path) {
  header <- "digraph nodes { \n node [shape=box] ;"
  footer <- "}"
  body <- create_dot_body(tree, include.na.path = include.na.path)

  dot_string <- paste(header, body, footer, sep = "\n")

  return(dot_string)
}

#' Plot a GRF tree object.
#'
#' The direction NAs are sent are indicated with the arrow fill. An empty arrow indicates
#' that NAs are sent that way. If trained without missing values, both arrows are filled.
#' @param x The tree to plot
#' @param include.na.path A boolean toggling whether to include the path of missing values or not.
#'  It defaults to whether the forest was trained with NAs.
#' @param ... Additional arguments (currently ignored).
#'
#' @method plot grf_tree
#' @examples
#' \dontrun{
#' # Plot a tree in the forest (requires the `DiagrammeR` package).
#' n <- 500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' c.forest <- causal_forest(X, Y, W)
#' plot(tree <- get_tree(c.forest, 1))
#' # Saving a plot in .svg can be done with the `DiagrammeRsvg` package.
#' install.packages("DiagrammeRsvg")
#' tree.plot = plot(tree)
#' cat(DiagrammeRsvg::export_svg(tree.plot), file = 'plot.svg')
#' }
#' @export
plot.grf_tree <- function(x, include.na.path = NULL, ...) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package \"DiagrammeR\" must be installed to plot trees.")
  }

  if (is.null(include.na.path)) {
    include.na.path <- x[["has.missing.values"]]
  }

  dot_file <- export_graphviz(x, include.na.path = include.na.path)
  DiagrammeR::grViz(dot_file)
}
