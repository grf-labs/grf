#' Print a generalized random forest.
print.grf <- function(forest, decay.exponent=0.5, max.depth=12) {
    split.freq = compute_split_frequencies(forest, max.depth)
    split.freq = split.freq / pmax(1, rowSums(split.freq))

    weight = (1:nrow(split.freq))^decay.exponent
    var.importance = t(split.freq) %*% weight / sum(weight)
    var.importance = c(round(var.importance, 3))
    names(var.importance) = 1:length(var.importance)

    main.class = class(forest)[1]
    num.observations = ncol(forest$original.data)

    print(paste("GRF object of type", main.class))
    print(paste("Number of observations: ", num.observations))
 
    print("Variable importance:")
    print(var.importance)
}
