#' Print a GRF forest object.
#' @param x The tree to print.
#' @param decay.exponent A tuning parameter that controls the importance of split depth.
#' @param max.depth The maximum depth of splits to consider.
#' @param ... Additional arguments (currently ignored).
#' @export
print.grf <- function(x, decay.exponent=2, max.depth=4, ...) {
    var.importance = variable_importance(x, decay.exponent, max.depth)
    var.importance = c(round(var.importance, 3))
    names(var.importance) = 1:length(var.importance)

    main.class = class(x)[1]
    num.samples= nrow(x$X.orig)

    cat("GRF forest object of type", main.class, "\n")
    cat("Number of trees: ", x$num.trees, "\n")
    cat("Number of training samples:", num.samples, "\n")
 
    cat("Variable importance:", "\n")
    print(var.importance)
}

#' Print a GRF tree object.
#' @param x The tree to print.
#' @param ... Additional arguments (currently ignored).
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

create_dot_body <- function(nodes, index=1) {
    
    node <- nodes[[index]]
    
    # Leaf case: print label only
    if (node$is_leaf) {
        num_samples = length(node$samples)
        line_label <- paste(index-1, ' [shape=box,style=filled,color=".7 .3 1.0" , label="leaf' , index,' \nsamples = ',num_samples,'"];')
        return(line_label)
    }
    
    # Non-leaf case: print label, child edges
    if(!is.null(node$left_child)){
        edge = paste(index-1, "->", node$left_child-1)
        if(index ==1 ){
            edge_info_left <- paste(edge, '[labeldistance=2.5, labelangle=45, headlabel="True"];')
        }
        else{
            edge_info_left <- paste(edge, " ;")
        }
    } 
    else{edge_info_right<-NULL}
    
    if(!is.null(node$right_child)){
        edge = paste(index-1, "->", node$right_child-1)
        if(index==1){
            edge_info_right <- paste(edge, '[labeldistance=2.5, labelangle=-45, headlabel="False"]')
        }
        else{
            edge_info_right <- paste(edge, " ;")
        }
    }else{edge_info_right<-NULL}
    
    node_info <- paste(index-1, '[label="split_variable', node$split_variable, '<=', round(node$split_value,2), '"] ;')
    
    this_lines <- paste(node_info,
                        edge_info_left,
                        edge_info_right, sep="\n")
    
    left_child_lines <- ifelse(!is.null(node$left_child),
                               create_dot_body(nodes, index=node$left_child),
                               NULL)
    
    right_child_lines <- ifelse(!is.null(node$right_child),
                                create_dot_body(nodes, index=node$right_child),
                                NULL)
    
    lines <- paste(this_lines, left_child_lines, right_child_lines, sep="\n")
    
    return(lines)
}

#' Plot a GRF tree object.
#' @param x The tree to print.
#' @param ... Additional arguments (currently ignored).
#' @export
plot.grf_tree <- function(x, ...){
    nodes=tree$nodes
    
    header <- "digraph nodes { \n node [shape=box] ;"
    footer <- "}"
    body <- create_dot_body(nodes)
    
    dot_string <- paste(header, body, footer, sep="\n")
    
    DiagrammeR::grViz(dot_string)
}
