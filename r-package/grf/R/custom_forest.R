#' Custom forest
#' 
#' Trains a custom forest model.
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty is used, these subsamples will
#'                        further be cut in half.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param min.node.size Minimum number of observations in each tree leaf.
#' @param keep.inbag Currently not used.
#' @param honesty Should honest splitting (i.e., sub-sample splitting) be used?   
#' @param alpha Maximum imbalance of a split.   
#' @param seed The seed of the c++ random number generator.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A trained regression forest object.
#' @export
custom_forest <- function(X, Y, sample.fraction = 0.5, mtry = ceiling(2*ncol(X)/3), 
    num.trees = 2000, num.threads = NULL, min.node.size = NULL, keep.inbag = FALSE, 
    honesty = TRUE, alpha = 0.05, seed = NULL) {
    
    sparse.data <- as.matrix(0)
    
    if (is.null(mtry)) {
        mtry <- 0
    } else if (!is.numeric(mtry) | mtry < 0) {
        stop("Error: Invalid value for mtry")
    }
    
    verbose = FALSE
    
    if (is.null(num.threads)) {
        num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
    }
    
    if (is.null(min.node.size)) {
        min.node.size <- 0
    } else if (!is.numeric(min.node.size) | min.node.size < 0) {
        stop("Error: Invalid value for min.node.size")
    }
    
    sample.with.replacement <- FALSE
    
    if (!is.logical(keep.inbag)) {
        stop("Error: Invalid value for keep.inbag")
    }
    
    if (!is.numeric(sample.fraction) | sample.fraction <= 0 | sample.fraction > 1) {
        stop("Error: Invalid value for sample.fraction. Please give a value in (0,1].")
    }
    
    if (is.null(seed)) {
        seed <- runif(1, 0, .Machine$integer.max)
    }
    
    input.data <- as.matrix(cbind(X, Y))
    variable.names <- c(colnames(X), "outcome")
    outcome.index <- ncol(input.data)
    no.split.variables <- numeric(0)
    ci.group.size <- 1
    
    forest <- custom_train(input.data, outcome.index, sparse.data,
        variable.names, mtry, num.trees, verbose, num.threads, min.node.size, sample.with.replacement,
        keep.inbag, sample.fraction, no.split.variables, seed, honesty, ci.group.size, alpha)
    
    forest[["original.data"]] <- input.data
    class(forest) <- c("custom_forest", "grf")
    forest
}

#' Predict with a custom forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions.
#' @export
predict.custom_forest <- function(object, newdata = NULL, num.threads = NULL, ...) {
    
    if (is.null(num.threads)) {
        num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
    }
    
    sparse.data <- as.matrix(0)
    variable.names <- character(0)
    
    forest.short <- object[-which(names(object) == "original.data")]
    
    if (!is.null(newdata)) {
        input.data <- as.matrix(cbind(newdata, NA))
        custom_predict(forest.short, input.data, sparse.data, variable.names, 
            num.threads)
    } else {
        input.data <- object[["original.data"]]
        custom_predict_oob(forest.short, input.data, sparse.data, variable.names, 
            num.threads)
    }
}
