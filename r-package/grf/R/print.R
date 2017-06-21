#' Print a GRF forest object.
print.grf <- function(forest, decay.exponent=0.5, max.depth=12) {
    split.freq = compute_split_frequencies(forest, max.depth)
    split.freq = split.freq / pmax(1, rowSums(split.freq))

    weight = (1:nrow(split.freq))^decay.exponent
    var.importance = t(split.freq) %*% weight / sum(weight)
    var.importance = c(round(var.importance, 3))
    names(var.importance) = 1:length(var.importance)

    main.class = class(forest)[1]
    num.samples= ncol(forest$original.data)

    cat("GRF forest object of type", main.class, "\n")
    cat("Number of training samples:", num.samples, "\n")
 
    cat("Variable importance:", "\n")
    print(var.importance)
}

#' Print a GRF tree object.
print.grf.tree <- function(tree) {
    cat("GRF tree object", "\n")
    cat("Number of training samples: ", tree$num_samples, "\n")
    cat("Variable splits:", "\n")

    # Add the index of each node as an attribute for easy access.
    nodes = lapply(1:length(tree$nodes), function(i) {
        node = tree$nodes[[i]]
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
            split.var.name = if (tree$columns[split.var] != "") tree$columns[split.var] else paste("X", split.var, sep=".")
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
