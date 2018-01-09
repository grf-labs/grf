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
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#' @param honesty Whether or not honest splitting (i.e., sub-sample splitting) should be used.
#' @param alpha Maximum imbalance of a split.   
#' @param seed The seed for the C++ random number generator.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A trained regression forest object.
#'
#' @examples \dontrun{
#' # Train a custom forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' c.forest = custom_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#' }
#'
#' @export
custom_forest <- function(X, Y, sample.fraction = 0.5, mtry = NULL, 
    num.trees = 2000, num.threads = NULL, min.node.size = NULL,
    honesty = TRUE, alpha = 0.05, seed = NULL) {

    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    
    mtry <- validate_mtry(mtry, X)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)
    
    no.split.variables <- numeric(0)
    sample.with.replacement <- FALSE
    verbose <- FALSE
    keep.inbag <- FALSE
    
    data <- create_data_matrices(X, Y)
    variable.names <- c(colnames(X), "outcome")
    outcome.index <- ncol(X) + 1
    no.split.variables <- numeric(0)
    ci.group.size <- 1
    
    forest <- custom_train(data$default, data$sparse, outcome.index,
        variable.names, mtry, num.trees, verbose, num.threads, min.node.size, sample.with.replacement,
        keep.inbag, sample.fraction, no.split.variables, seed, honesty, ci.group.size, alpha)
    
    forest[["X.orig"]] <- X
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
#'
#' @examples \dontrun{
#' # Train a custom forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' c.forest = custom_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#' }
#'
#' @export
predict.custom_forest <- function(object, newdata = NULL, num.threads = NULL, ...) {
    
    if (is.null(num.threads)) {
        num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
    }
    
    variable.names <- character(0)
    
    forest.short <- object[-which(names(object) == "X.orig")]
    
    if (!is.null(newdata)) {
        data <- create_data_matrices(newdata, NA)
        custom_predict(forest.short, data$default, data$sparse,
                       variable.names, num.threads)
    } else {
        data <- create_data_matrices(object[["X.orig"]], NA)
        custom_predict_oob(forest.short, data$default, data$sparse,
                           variable.names, num.threads)
    }
}
