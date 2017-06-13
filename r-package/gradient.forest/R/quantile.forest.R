#' Quantile forest
#' 
#' Trains a regression forest that can be used to estimate
#' quantiles of the conditional distribution of Y given X = x.
#'
#' @param X The covariates used in the quantile regression.
#' @param Y The outcome.
#' @param quantiles Vector of quantiles used to calibrate the forest.
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
#' @param ci.group.size The forst will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. [Note: confidence intervals for quantile forests
#'                      are not yet implemented.]
#' @param seed The seed of the c++ random number generator.
#'
#' @return A trained quantile forest object.
#' @export
quantile.forest <- function(X, Y, quantiles = c(0.1, 0.5, 0.9), sample.fraction = 0.5, 
    mtry = ceiling(2*ncol(X)/3), num.trees = 2000, num.threads = NULL, min.node.size = NULL, 
    keep.inbag = FALSE, seed = NULL, ci.group.size = 2, alpha = 0.10, honesty = TRUE) {
    
    if (!is.numeric(quantiles) | length(quantiles) < 1) {
        stop("Error: Must provide numeric quantiles")
    } else if (min(quantiles) <= 0 | max(quantiles) >= 1) {
        stop("Error: Quantiles must be in (0, 1)")
    }
    
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
    outcome.index.zeroindexed <- outcome.index - 1
    no.split.variables <- numeric(0)
    
    forest <- quantile_train(quantiles, input.data, outcome.index.zeroindexed, sparse.data, 
        variable.names, mtry, num.trees, verbose, num.threads, min.node.size, sample.with.replacement, 
        keep.inbag, sample.fraction, no.split.variables, seed, honesty, ci.group.size, alpha)
    
    forest[["original.data"]] <- input.data
    class(forest) <- c("quantile.forest", "grf")
    forest
}

#' Predict with a quantile forest
#' 
#' Gets estimates of the conditional quantiles of Y given X using a trained forest.
#'
#' @param forest The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param quantiles Vector of quantiles at which estimates are required.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#'
#' @return Predictions for each test point and each desired quantile.
#' @export

predict.quantile.forest <- function(forest, newdata = NULL, quantiles = c(0.1, 0.5, 
    0.9), num.threads = NULL) {
    
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
    
    forest.short <- forest[-which(names(forest) == "original.data")]
    
    if (!is.null(newdata)) {
        input.data <- as.matrix(cbind(newdata, NA))
        quantile_predict(forest, quantiles, input.data, sparse.data, variable.names, 
            num.threads)
    } else {
        input.data <- forest[["original.data"]]
        quantile_predict_oob(forest, quantiles, input.data, sparse.data, variable.names, 
            num.threads)
    }
}
