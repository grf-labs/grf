#' Quantile forest
#' 
#' Trains a regression forest that can be used to estimate
#' quantiles of the conditional distribution of Y given X = x.
#'
#' @param X The covariates used in the quantile regression.
#' @param Y The outcome.
#' @param quantiles Vector of quantiles used to calibrate the forest.
#' @param regression.splitting Whether to use regression splits when growing trees instead
#'                             of specialized splits based on the quantiles (the default).
#'                             Setting this flag to true corresponds to the approach to
#'                             quantile forests from Meinshausen (2006).
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
#' @param seed The seed of the c++ random number generator.
#' @param honesty Should honest splitting (i.e., sub-sample splitting) be used?
#' @param alpha Maximum imbalance of a split.
#'
#' @return A trained quantile forest object.
#' @export
quantile_forest <- function(X, Y, quantiles = c(0.1, 0.5, 0.9), regression.splitting = FALSE,
                            sample.fraction = 0.5, mtry = ceiling(2*ncol(X)/3), num.trees = 2000,
                            num.threads = NULL, min.node.size = NULL, seed = NULL, alpha = 0.05,
                            honesty = TRUE) {
    
    if (!is.numeric(quantiles) | length(quantiles) < 1) {
        stop("Error: Must provide numeric quantiles")
    } else if (min(quantiles) <= 0 | max(quantiles) >= 1) {
        stop("Error: Quantiles must be in (0, 1)")
    }
    
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    
    mtry <- validate_mtry(mtry)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)
    
    sparse.data <- as.matrix(0)
    no.split.variables <- numeric(0)
    sample.with.replacement <- FALSE
    verbose <- FALSE
    keep.inbag <- FALSE
    
    input.data <- as.matrix(cbind(X, Y))
    variable.names <- c(colnames(X), "outcome")
    outcome.index <- ncol(input.data)

    ci.group.size <- 1
    
    forest <- quantile_train(quantiles, regression.splitting, input.data, outcome.index, sparse.data,
        variable.names, mtry, num.trees, verbose, num.threads, min.node.size, sample.with.replacement,
        keep.inbag, sample.fraction, no.split.variables, seed, honesty, ci.group.size, alpha)
    
    forest[["original.data"]] <- input.data
    forest[["feature.indices"]] <- 1:ncol(X)
    class(forest) <- c("quantile_forest", "grf")
    forest
}

#' Predict with a quantile forest
#' 
#' Gets estimates of the conditional quantiles of Y given X using a trained forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param quantiles Vector of quantiles at which estimates are required.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Predictions for each test point and each desired quantile.
#' @export

predict.quantile_forest <- function(object,
                                    newdata = NULL,
                                    quantiles = c(0.1, 0.5, 0.9),
                                    num.threads = NULL, ...) {
    
    if (!is.numeric(quantiles) | length(quantiles) < 1) {
        stop("Error: Must provide numeric quantiles")
    } else if (min(quantiles) <= 0 | max(quantiles) >= 1) {
        stop("Error: Quantiles must be in (0, 1)")
    }
    
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
        quantile_predict(forest.short, quantiles, input.data, sparse.data, variable.names, 
                         num.threads)
    } else {
        input.data <- object[["original.data"]]
        quantile_predict_oob(forest.short, quantiles, input.data, sparse.data, variable.names, 
                             num.threads)
    }
}
