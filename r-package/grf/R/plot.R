#' Writes each node information 
#' If it is a leaf node: show it in different color, show number of samples, show leaf id
#' If it is a non-leaf node: show its splitting variable and splitting value
#' @param tree the tree to convert
#' @param index the index of the current node
create_dot_body <- function(tree, index=1) {
  
  node <- tree$nodes[[index]]
  
  # Leaf case: print label only
  if (node$is_leaf) {
    num_samples = length(node$samples)
    line_label <- paste(index-1, ' [shape=box,style=filled,color=".7 .3 1.0" , label="leaf node', ' \nsize = ',num_samples,'"];')
    return(line_label)
  }
  
  # Non-leaf case: print label, child edges
  if (!is.null(node$left_child)) {
    edge = paste(index-1, "->", node$left_child-1)
    if (index == 1){
      edge_info_left <- paste(edge, '[labeldistance=2.5, labelangle=45, headlabel="True"];')
    }
    else {
      edge_info_left <- paste(edge, " ;")
    }
  } 
  else{edge_info_right<-NULL}
  
  if (!is.null(node$right_child)) {
    edge = paste(index-1, "->", node$right_child-1)
    if (index==1) {
      edge_info_right <- paste(edge, '[labeldistance=2.5, labelangle=-45, headlabel="False"]')
    } else {
      edge_info_right <- paste(edge, " ;")
    }
  } else {
    edge_info_right<-NULL
  }
  
  variable_name = tree$columns[node$split_variable]
  node_info <- paste(index-1, '[label="', variable_name, '<=', round(node$split_value,2), '"] ;')
  
  this_lines <- paste(node_info,
                      edge_info_left,
                      edge_info_right, sep="\n")
  
  left_child_lines <- ifelse(!is.null(node$left_child),
                             create_dot_body(tree, index=node$left_child),
                             NULL)
  
  right_child_lines <- ifelse(!is.null(node$right_child),
                              create_dot_body(tree, index=node$right_child),
                              NULL)
  
  lines <- paste(this_lines, left_child_lines, right_child_lines, sep="\n")
  
  return(lines)
}

#' Export a tree in DOT format.
#' This function generates a GraphViz representation of the tree,
#' which is then written into `dot_string`. 
#' @param tree the tree to convert
export_graphviz <- function(tree){
  header <- "digraph nodes { \n node [shape=box] ;"
  footer <- "}"
  body <- create_dot_body(tree)
  
  dot_string <- paste(header, body, footer, sep="\n")
  
  return(dot_string)
}

#' Plot a GRF tree object.
#' @param x The tree to plot
#' @param ... Additional arguments (currently ignored).
#'
#' @method plot grf_tree
#' @export
plot.grf_tree <- function(x, ...){
  dot_file <- export_graphviz(x)
  DiagrammeR::grViz(dot_file)
}

