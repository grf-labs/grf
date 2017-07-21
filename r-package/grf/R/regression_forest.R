#' Regression forest
#' 
#' Trains a regression forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]
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
#' @param honesty Should honest splitting (i.e., sub-sample splitting) be used?
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2.
#' @param alpha Maximum imbalance of a split.
#' @param lambda A tuning parameter to control the amount of split regularization (experimental).
#' @param downweight.penalty Whether or not the regularization penalty should be downweighted (experimental).
#' @param seed The seed of the c++ random number generator.
#'
#' @return A trained regression forest object.
#'
#' @examples
#' # Train a standard regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.forest = regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' r.pred = predict(r.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred = predict(r.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' r.forest = regression_forest(X, Y, num.trees = 100)
#' r.pred = predict(r.forest, X.test, estimate.variance = TRUE)
#'
#' @export
regression_forest <- function(X, Y, sample.fraction = 0.5, mtry = ceiling(2*ncol(X)/3), 
                              num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                              honesty = TRUE, ci.group.size = 2, alpha = 0.05, lambda = 0.0,
                              downweight.penalty = FALSE, seed = NULL) {
    
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    
    mtry <- validate_mtry(mtry)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)
    
    no.split.variables <- numeric(0)
    sample.with.replacement <- FALSE
    verbose <- FALSE
    keep.inbag <- FALSE
    
    input.data <- as.matrix(cbind(X, Y))
    variable.names <- c(colnames(X), "outcome")
    outcome.index <- ncol(input.data)
    
    forest <- regression_train(input.data, outcome.index, variable.names, mtry, num.trees,
        verbose, num.threads, min.node.size, sample.with.replacement, keep.inbag, sample.fraction,
        no.split.variables, seed, honesty, ci.group.size, alpha, lambda, downweight.penalty)
    
    forest[["ci.group.size"]] <- ci.group.size
    forest[["original.data"]] <- input.data
    forest[["feature.indices"]] <- 1:ncol(X)
    class(forest) <- c("regression_forest", "grf")
    forest
}

#' Predict with a regression forest
#' 
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A vector of predictions.
#'
#' @examples
#' # Train a standard regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.forest = regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' r.pred = predict(r.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred = predict(r.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' r.forest = regression_forest(X, Y, num.trees = 100)
#' r.pred = predict(r.forest, X.test, estimate.variance = TRUE)
#'
#' @export
predict.regression_forest <- function(object, newdata = NULL,
                                      num.threads = NULL,
                                      estimate.variance = FALSE,
                                      ...) {
    
    if (is.null(num.threads)) {
        num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
    }
    
    variable.names <- character(0)
    
    if (estimate.variance) {
        ci.group.size = object$ci.group.size
    } else {
        ci.group.size = 1
    }
    
    forest.short <- object[-which(names(object) == "original.data")]
    
    if (!is.null(newdata)) {
        input.data <- as.matrix(cbind(newdata, NA))
        regression_predict(forest.short, input.data, variable.names, 
                           num.threads, ci.group.size)
    } else {
        input.data <- object[["original.data"]]
        regression_predict_oob(forest.short, input.data, variable.names, 
                               num.threads, ci.group.size)
    }
}
