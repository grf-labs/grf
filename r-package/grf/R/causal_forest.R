#' Causal forest
#' 
#' Trains a causal forest that can be used to estimate
#' conditional average treatment effects tau(X). When
#' the treatment assignmnet W is binary and unconfounded,
#' we have tau(X) = E[Y(1) - Y(0) | X = x], where Y(0) and
#' Y(1) are potential outcomes corresponding to the two possible
#' treatment states. When W is continuous, we effectively estimate
#' an average partical effect Cov[Y, W | X = x] / Var[W | X = x],
#' and interpret it as a treatment effect given unconfoundedness.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
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
#' @param ci.group.size The forst will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2.
#' @param precompute.nuisance Should we first run regression forests to estimate
#'                            y(x) = E[Y|X=x] and w(x) = E[W|X=x], and then run a
#'                            causal forest on the residuals? This approach is
#'                            recommended, computational resources permitting.
#' @param alpha Maximum imbalance of a split.
#' @param seed The seed of the c++ random number generator.
#'
#' @return A trained causal forest object.
#' @export
causal_forest <- function(X, Y, W, sample.fraction = 0.5, mtry = ceiling(2*ncol(X)/3), 
                          num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                          honesty = TRUE, ci.group.size = 2, precompute.nuisance = TRUE,
                          alpha = 0.05, seed = NULL) {
    
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    if(length(W) != nrow(X)) { stop("W has incorrect length.") }
    
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
    
    split.regularization <- 0
    
    if (!precompute.nuisance) {
        
        input.data <- as.matrix(cbind(X, Y, W))
        Y.hat <- NULL
        W.hat <- NULL
        
    } else {
        
        forest.Y <- regression_forest(X, Y, sample.fraction = sample.fraction, mtry = mtry, 
                                      num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, 
                                      honesty = TRUE, seed = seed, ci.group.size = 1, alpha = alpha)
        Y.hat <- predict(forest.Y)$predictions
        
        forest.W <- regression_forest(X, W, sample.fraction = sample.fraction, mtry = mtry, 
                                      num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, 
                                      honesty = TRUE, seed = seed, ci.group.size = 1, alpha = alpha)
        W.hat <- predict(forest.W)$predictions
        
        input.data <- as.matrix(cbind(X, Y - Y.hat, W - W.hat))
        
    }
    
    
    variable.names <- c(colnames(X), "outcome", "treatment")
    outcome.index <- ncol(X) + 1
    treatment.index <- ncol(X) + 2
    instrument.index <- treatment.index
    
    forest <- instrumental_train(input.data, outcome.index, treatment.index,
        instrument.index, sparse.data, variable.names, mtry, num.trees, verbose,
        num.threads, min.node.size, sample.with.replacement, keep.inbag, sample.fraction,
        no.split.variables, seed, honesty, ci.group.size, split.regularization, alpha)
    
    forest[["ci.group.size"]] <- ci.group.size
    forest[["original.data"]] <- input.data
    forest[["feature.indices"]] <- 1:ncol(X)
    forest[["Y.orig"]] <- Y
    forest[["W.orig"]] <- W
    forest[["Y.hat"]] <- Y.hat
    forest[["W.hat"]] <- W.hat
    
    class(forest) <- c("causal_forest", "grf")
    forest
}

#' Predict with a causal forest
#' 
#' Gets estimates of tau(x) using a trained causal forest.
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
#' @return Vector of predictions, along with (optional) variance estimates.
#' @export
predict.causal_forest <- function(object, newdata = NULL, num.threads = NULL, estimate.variance = FALSE, ...) {
    predict.instrumental_forest(object, newdata, num.threads, estimate.variance)
}
