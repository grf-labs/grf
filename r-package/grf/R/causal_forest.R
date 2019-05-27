#' Causal forest
#'
#' Trains a causal forest that can be used to estimate
#' conditional average treatment effects tau(X). When
#' the treatment assignment W is binary and unconfounded,
#' we have tau(X) = E[Y(1) - Y(0) | X = x], where Y(0) and
#' Y(1) are potential outcomes corresponding to the two possible
#' treatment states. When W is continuous, we effectively estimate
#' an average partial effect Cov[Y, W | X = x] / Var[W | X = x],
#' and interpret it as a treatment effect given unconfoundedness.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome (must be a numeric vector with no NAs).
#' @param W The treatment assignment (must be a binary or real numeric vector with no NAs).
#' @param Y.hat Estimates of the expected responses E[Y | Xi], marginalizing
#'              over treatment. If Y.hat = NULL, these are estimated using
#'              a separate regression forest. See section 6.1.1 of the GRF paper for
#'              further discussion of this quantity.
#' @param W.hat Estimates of the treatment propensities E[W | Xi]. If W.hat = NULL,
#'              these are estimated using a separate regression forest.
#' @param sample.weights (experimental) Weights given to each sample in estimation.
#'                       If NULL, each observation receives the same weight.
#'                       Note: To avoid introducing confounding, weights should be
#'                       independent of the potential outcomes given X.
#' @param orthog.boosting (experimental) If TRUE, then when Y.hat = NULL or W.hat is NULL,
#'                 the missing quantities are estimated using boosted regression forests.
#'                 The number of boosting steps is selected automatically.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting).
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. When using the defaults (honesty = TRUE and
#'                         honesty.fraction = NULL), half of the data will be used for determining splits
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized.
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#'                         determining the imbalance of a split.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#' @param samples.per.cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster when training a tree. If NULL, we set samples.per.cluster to the size
#'                            of the smallest cluster. If some clusters are smaller than samples.per.cluster,
#'                            the whole cluster is used every time the cluster is drawn. Note that
#'                            clusters with less than samples.per.cluster observations get relatively
#'                            smaller weight than others in training the forest, i.e., the contribution
#'                            of a given cluster to the final forest scales with the minimum of
#'                            the number of observations in the cluster and samples.per.cluster.
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if false
#'                        NULL parameters are set to defaults.
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param num.fit.reps The number of forests used to fit the tuning model.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters.
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained causal forest object. If tune.parameters is enabled,
#'  then tuning information will be included through the `tuning.output` attribute.
#'
#' @examples \dontrun{
#' # Train a causal forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' c.forest = causal_forest(X, Y, W)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' c.pred = predict(c.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' c.forest = causal_forest(X, Y, W, num.trees = 4000)
#' c.pred = predict(c.forest, X.test, estimate.variance = TRUE)
#'
#' # In some examples, pre-fitting models for Y and W separately may
#' # be helpful (e.g., if different models use different covariates).
#' # In some applications, one may even want to get Y.hat and W.hat
#' # using a completely different method (e.g., boosting).
#' n = 2000; p = 20
#' X = matrix(rnorm(n * p), n, p)
#' TAU = 1 / (1 + exp(-X[, 3]))
#' W = rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
#' Y = pmax(X[, 2] + X[, 3], 0) + rowMeans(X[, 4:6]) / 2 + W * TAU + rnorm(n)
#'
#' forest.W = regression_forest(X, W, tune.parameters = TRUE)
#' W.hat = predict(forest.W)$predictions
#'
#' forest.Y = regression_forest(X, Y, tune.parameters = TRUE)
#' Y.hat = predict(forest.Y)$predictions
#'
#' forest.Y.varimp = variable_importance(forest.Y)
#'
#' # Note: Forests may have a hard time when trained on very few variables
#' # (e.g., ncol(X) = 1, 2, or 3). We recommend not being too aggressive
#' # in selection.
#' selected.vars = which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)
#'
#' tau.forest = causal_forest(X[,selected.vars], Y, W,
#'                            W.hat = W.hat, Y.hat = Y.hat,
#'                            tune.parameters = TRUE)
#' tau.hat = predict(tau.forest)$predictions
#' }
#'
#' @export
causal_forest <- function(X, Y, W,
                          Y.hat = NULL,
                          W.hat = NULL,
                          sample.weights = NULL,
                          orthog.boosting = FALSE,
                          sample.fraction = 0.5,
                          mtry = NULL,
                          num.trees = 2000,
                          min.node.size = NULL,
                          honesty = TRUE,
                          honesty.fraction = NULL,
                          ci.group.size = 2,
                          alpha = NULL,
                          imbalance.penalty = NULL,
                          stabilize.splits = TRUE,
                          clusters = NULL,
                          samples.per.cluster = NULL,
                          tune.parameters = FALSE,
                          num.fit.trees = 200,
                          num.fit.reps = 50,
                          num.optimize.reps = 1000,
                          compute.oob.predictions = TRUE,
                          num.threads = NULL,
                          seed = NULL) {
    validate_X(X)
    validate_sample_weights(sample.weights, X)
    Y = validate_observations(Y, X)
    W = validate_observations(W, X)

    num.threads <- validate_num_threads(num.threads)
    seed <- validate_seed(seed)
    clusters <- validate_clusters(clusters, X)
    samples.per.cluster <- validate_samples_per_cluster(samples.per.cluster, clusters)
    honesty.fraction <- validate_honesty_fraction(honesty.fraction, honesty)

    reduced.form.weight <- 0

    if (is.null(Y.hat) && !orthog.boosting) {
      forest.Y <- regression_forest(X, Y, sample.weights = sample.weights, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples.per.cluster = samples.per.cluster);
      Y.hat <- predict(forest.Y)$predictions
    } else if (is.null(Y.hat) && orthog.boosting) {
      forest.Y <- boosted_regression_forest(X, Y, sample.weights = sample.weights, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples.per.cluster = samples.per.cluster);
      Y.hat <- predict(forest.Y)$predictions
    } else if (length(Y.hat) == 1) {
      Y.hat <- rep(Y.hat, nrow(X))
    } else if (length(Y.hat) != nrow(X)) {
      stop("Y.hat has incorrect length.")
    }

    if (is.null(W.hat) && !orthog.boosting) {
      forest.W <- regression_forest(X, W, sample.weights = sample.weights, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples.per.cluster = samples.per.cluster);
      W.hat <- predict(forest.W)$predictions

    } else if (is.null(W.hat) && orthog.boosting) {
      forest.W <- boosted_regression_forest(X, W, sample.weights = sample.weights, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples.per.cluster = samples.per.cluster);
      W.hat <- predict(forest.W)$predictions
    } else if (length(W.hat) == 1) {
      W.hat <- rep(W.hat, nrow(X))
    } else if (length(W.hat) != nrow(X)) {
      stop("W.hat has incorrect length.")
    }

    if (tune.parameters) {
      tuning.output <- tune_causal_forest(X, Y, W, Y.hat, W.hat,
                                          sample.weights = sample.weights,
                                          num.fit.trees = num.fit.trees,
                                          num.fit.reps = num.fit.reps,
                                          num.optimize.reps = num.optimize.reps,
                                          min.node.size = min.node.size,
                                          sample.fraction = sample.fraction,
                                          mtry = mtry,
                                          alpha = alpha,
                                          imbalance.penalty = imbalance.penalty,
                                          stabilize.splits = stabilize.splits,
                                          num.threads = num.threads,
                                          honesty = honesty,
                                          honesty.fraction = honesty.fraction,
                                          seed = seed,
                                          clusters = clusters,
                                          samples.per.cluster = samples.per.cluster)
      tunable.params <- tuning.output$params
    } else {
      tunable.params <- c(
        min.node.size = validate_min_node_size(min.node.size),
        sample.fraction = validate_sample_fraction(sample.fraction),
        mtry = validate_mtry(mtry, X),
        alpha = validate_alpha(alpha),
        imbalance.penalty = validate_imbalance_penalty(imbalance.penalty))
    }

    Y.centered = Y - Y.hat
    W.centered = W - W.hat

    data <- create_data_matrices(X, Y.centered, W.centered, sample.weights = sample.weights)
    outcome.index <- ncol(X) + 1
    treatment.index <- ncol(X) + 2
    sample.weight.index <- ncol(X) + 3

    forest <- causal_train(data$default, data$sparse,
                           outcome.index, treatment.index, sample.weight.index,
                           !is.null(sample.weights),
                           as.numeric(tunable.params["mtry"]),
                           num.trees,
                           as.numeric(tunable.params["min.node.size"]),
                           as.numeric(tunable.params["sample.fraction"]),
                           honesty,
                           coerce_honesty_fraction(honesty.fraction),
                           ci.group.size,
                           reduced.form.weight,
                           as.numeric(tunable.params["alpha"]),
                           as.numeric(tunable.params["imbalance.penalty"]),
                           stabilize.splits,
                           clusters,
                           samples.per.cluster,
                           compute.oob.predictions,
                           num.threads,
                           seed)

    class(forest) <- c("causal_forest", "grf")
    forest[["ci.group.size"]] <- ci.group.size
    forest[["X.orig"]] <- X
    forest[["Y.orig"]] <- Y
    forest[["W.orig"]] <- W
    forest[["Y.hat"]] <- Y.hat
    forest[["W.hat"]] <- W.hat
    forest[["clusters"]] <- clusters
    forest[["tunable.params"]] <- tunable.params
    forest[["sample.weights"]] <- sample.weights
    if (tune.parameters)
      forest[["tuning.output"]] <- tuning.output
    forest
}

#' Predict with a causal forest
#'
#' Gets estimates of tau(x) using a trained causal forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param linear.correction.variables Optional subset of indexes for variables to be used in local
#'                   linear prediction. If NULL, standard GRF prediction is used. Otherwise,
#'                   we run a locally weighted linear regression on the included variables.
#'                   Please note that this is a beta feature still in development, and may slow down
#'                   prediction considerably. Defaults to NULL.
#' @param ll.lambda Ridge penalty for local linear predictions
#' @param ll.weight.penalty Option to standardize ridge penalty by covariance (TRUE),
#'                  or penalize all covariates equally (FALSE). Penalizes equally by default.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions, along with estimates of the error and
#'         (optionally) its variance estimates. Column 'predictions' contains estimates
#'         of the conditional average treatent effect (CATE). The square-root of
#'         column 'variance.estimates' is the standard error of CATE.
#'         For out-of-bag estimates, we also output the following error measures.
#'         First, column 'debiased.error' contains estimates of the 'R-loss' criterion,
#          a quantity that is related to the true (infeasible) mean-squared error
#'         (See Nie and Wager 2017 for a justification). Second, column 'excess.error'
#'         contains jackknife estimates of the Monte-carlo error (Wager, Hastie, Efron 2014),
#'         a measure of how unstable estimates are if we grow forests of the same size
#'         on the same data set. The sum of 'debiased.error' and 'excess.error' is the raw error
#'         attained by the current forest, and 'debiased.error' alone is an estimate of the error
#'         attained by a forest with an infinite number of trees. We recommend that users grow
#'         enough forests to make the 'excess.error' negligible.
#'
#' @examples \dontrun{
#' # Train a causal forest.
#' n = 100; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' c.forest = causal_forest(X, Y, W)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' c.pred = predict(c.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' c.forest = causal_forest(X, Y, W, num.trees = 500)
#' c.pred = predict(c.forest, X.test, estimate.variance = TRUE)
#' }
#'
#' @method predict causal_forest
#' @export
predict.causal_forest <- function(object, newdata = NULL,
                                  linear.correction.variables = NULL,
                                  ll.lambda = 0.1,
                                  ll.weight.penalty = FALSE,
                                  num.threads = NULL, estimate.variance = FALSE, ...) {

    # If possible, use pre-computed predictions.
    if (is.null(newdata) & !estimate.variance & !is.null(object$predictions) & is.null(linear.correction.variables)) {
        return(data.frame(predictions=object$predictions,
                          debiased.error=object$debiased.error,
                          excess.error=object$excess.error))
    }

    forest.short <- object[-which(names(object) == "X.orig")]

    X = object[["X.orig"]]
    Y.centered = object[["Y.orig"]] - object[["Y.hat"]]
    W.centered = object[["W.orig"]] - object[["W.hat"]]
    train.data <- create_data_matrices(X, Y.centered, W.centered)

    outcome.index <- ncol(X) + 1
    treatment.index <- ncol(X) + 2

    num.threads <- validate_num_threads(num.threads)

    local.linear = !is.null(linear.correction.variables)
    if(local.linear){
        linear.correction.variables = validate_ll_vars(linear.correction.variables, ncol(X))
        ll.lambda = validate_ll_lambda(ll.lambda)

        # subtract 1 to account for C++ indexing
        linear.correction.variables <- linear.correction.variables - 1
    }

    if (!is.null(newdata) ) {
        validate_newdata(newdata, object$X.orig)
        data = create_data_matrices(newdata)
        if (!local.linear) {
            ret <- causal_predict(forest.short, train.data$default, train.data$sparse,
                    outcome.index, treatment.index, data$default, data$sparse, num.threads, estimate.variance)
        } else {
            ret <- ll_causal_predict(forest.short, data$default, train.data$default, data$sparse, train.data$sparse,
                    outcome.index, treatment.index, ll.lambda, ll.weight.penalty, linear.correction.variables, num.threads)
        }
    } else {
        if (!local.linear) {
            ret <- causal_predict_oob(forest.short, train.data$default, train.data$sparse,
                    outcome.index, treatment.index, num.threads, estimate.variance)
        } else {
            ret <- ll_causal_predict_oob(forest.short, train.data$default, train.data$sparse,
                    outcome.index, treatment.index, ll.lambda, ll.weight.penalty, linear.correction.variables, num.threads)
        }
    }

    # Convert list to data frame.
    empty = sapply(ret, function(elem) length(elem) == 0)
    do.call(cbind.data.frame, ret[!empty])
}
