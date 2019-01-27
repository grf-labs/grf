#' Boosted regression forest
#'
#' Trains a boosted regression forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]. Selects
#' number of boosting iterations based on cross-validation.
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
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
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed.
#' @param seed The seed for the C++ random number generator.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#' @param samples_per_cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster when training a tree. If NULL, we set samples_per_cluster to the size
#'                            of the smallest cluster. If some clusters are smaller than samples_per_cluster,
#'                            the whole cluster is used every time the cluster is drawn. Note that
#'                            clusters with less than samples_per_cluster observations get relatively
#'                            smaller weight than others in training the forest, i.e., the contribution
#'                            of a given cluster to the final forest scales with the minimum of
#'                            the number of observations in the cluster and samples_per_cluster.
#' @param tune.parameters If true, NULL parameters are tuned by cross-validation; if false
#'                        NULL parameters are set to defaults.
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param num.fit.reps The number of forests used to fit the tuning model.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters.
#' @param num.steps The number of boosting iterations. If NULL, selected by cross-validation
#' @param max.steps The maximum number of boosting iterations
#'
#' @return A trained regression forest object.
#'
#' @examples \dontrun{
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
#' }
#'
#' @export
boosted_regression_forest <- function(X, Y,
                              sample.fraction = 0.5,
                              mtry = NULL,
                              num.trees = 2000,
                              num.threads = NULL,
                              min.node.size = NULL,
                              honesty = TRUE,
                              honesty.fraction = NULL,
                              ci.group.size = 2,
                              alpha = NULL,
                              imbalance.penalty = NULL,
                              compute.oob.predictions = TRUE,
                              seed = NULL,
                              clusters = NULL,
                              samples_per_cluster = NULL,
                              tune.parameters = FALSE,
                              num.fit.trees = 10,
                              num.fit.reps = 100,
                              num.optimize.reps = 1000,
                              num.steps = NULL,
                              max.steps = 4) {

    #TODO: Replace this with tuning procedure to find the right number of steps
    if(is.null(num.steps)) {
        num.steps = max.steps
    }
    boosted.forest = NULL
    boosted.forest[["forests"]] = list()
    forest.Y <- regression_forest(X, Y, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                  num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                  honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                clusters = clusters, samples_per_cluster = samples_per_cluster);
    Y.hat <- predict(forest.Y)$predictions
    boosted.forest$forests[[1]] <- forest.Y

    for (step in 2:max.steps) {
        E.hat <- Y - Y.hat
        forest.E <- regression_forest(X,E.hat, sample.fraction = sample.fraction, mtry = mtry, tune.parameters = tune.parameters,
                                      num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                      honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
        Y.hat <- Y.hat + predict(forest.E)$predictions
        boosted.forest[["forests"]][[step]] <- forest.E
    }

    if (compute.oob.predictions) {
        boosted.forest[["predictions"]] <- Y.hat
        #boosted.forest[["debiased.error"]] <- oob.pred$debiased.error
    }

    class(boosted.forest) <- c("boosted_regression_forest")
    boosted.forest
}

#' Predict with a boosted regression forest. No option for confidence intervals.
#'
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @return A vector of predictions.
#'
#' @examples \dontrun{
#' # Train a standard regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.forest = boosted_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' r.pred = predict(r.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred = predict(r.forest)
#'
#' }
#'
#' @method predict boosted_regression_forest
#' @export
predict.boosted_regression_forest <- function(object,
                                        newdata=NULL,
                                      num.threads = NULL
                                      ) {


    #confidence interval estimation not possible with boosted forests
    estimate.variance <- FALSE
    # If possible, use pre-computed predictions.
    if (is.null(newdata) & !is.null(object$predictions)) {
        return(data.frame(predictions=object$predictions))
        ##    ,debiased.error=object$debiased.error))
    }
    forests <- object[["forests"]]
    if (!is.null(newdata)) {
        ## TODO: Implement newdata option
    }
    else {
        Y <- forests[[1]][["Y.orig"]]
        Y.hat <- predict(forests[[1]])$predictions
        for (f in 2:length(forests)) {
            Y.hat <- Y.hat + predict(forests[[f]])$predictions
        }
    }
    data.frame(predictions=Y.hat)
}
