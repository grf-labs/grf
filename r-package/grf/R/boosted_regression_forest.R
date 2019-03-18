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
#' @param max.steps The maximum number of boosting iterations to try when num.steps is not specified
#' @param num.trees.tune If num.steps is NULL, the number of trees used to test a new boosting step when tuning
#'        num.steps
#'
#' @return A trained regression forest object.
#'
#' @examples \dontrun{
#' # Train a boosted regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.b.forest = boosted_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' r.pred = predict(r.b.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred = predict(r.b.forest)
#'
#' #Check how many boosting iterations were used
#' print(length(r.b.forest$forests))
#' }
#'
#' @export
boosted_regression_forest <- function(X, Y,
                              sample.fraction = 0.5,
                              mtry = NULL,
                              num.trees = 500,
                              num.threads = NULL,
                              min.node.size = NULL,
                              honesty = TRUE,
                              honesty.fraction = NULL,
                              ci.group.size = 2,
                              alpha = NULL,
                              imbalance.penalty = NULL,
                              seed = NULL,
                              clusters = NULL,
                              samples_per_cluster = NULL,
                              tune.parameters = FALSE,
                              num.fit.trees = 10,
                              num.fit.reps = 100,
                              num.optimize.reps = 1000,
                              num.steps = NULL,
                              max.steps = 6,
                              num.trees.tune = 10) {

    boosted.forest = NULL
    boosted.forest[["forests"]] = list()
    forest.Y <- regression_forest(X, Y, sample.fraction = sample.fraction,
                                mtry = mtry, tune.parameters = tune.parameters,
                                  num.trees = num.trees,
                                  num.threads = num.threads,
                                  min.node.size = min.node.size, honesty = honesty,
                                  honesty.fraction = honesty.fraction,
                                  seed = seed, ci.group.size = ci.group.size,
                                  alpha = alpha,
                                  imbalance.penalty = imbalance.penalty,
                                clusters = clusters, samples_per_cluster = samples_per_cluster);
    forest.Y.p <- predict(forest.Y)
    #save tuned parameters for use on future boosting iterations 
    tunable.params <- forest.Y$tunable.params
    Y.hat <- forest.Y.p$prediction
    error.debiased <- forest.Y.p$debiased.error
    boosted.forest$forests[[1]] <- forest.Y

    if(!is.null(num.steps)) {
        for (step in 2:num.steps) {
            E.hat <- Y - Y.hat
            forest.E <- regression_forest(X,E.hat, sample.fraction = tunable.params["sample.fraction"],
                                        mtry = tunable.params["mtry"], tune.parameters = FALSE,
                                          num.trees = num.trees,
                                           num.threads = num.threads,
                                           min.node.size = tunable.params["min.node.size"],
                                           honesty = honesty,
                                          honesty.fraction = honesty.fraction,
                                           seed = seed, ci.group.size = ci.group.size,
                                          alpha = tunable.params["alpha"],
                                          imbalance.penalty = tunable.params["imbalance.penalty"],
                                        clusters = clusters, samples_per_cluster = samples_per_cluster);
            current.pred <- predict(forest.E)
            Y.hat <- Y.hat + current.pred$predictions
            error.debiased <- current.pred$debiased.error
            boosted.forest[["forests"]][[step]] <- forest.E
        }
    }
    else {
        step <- 2
        tolerance <- 0.99
        still_boosting <- TRUE
        while(still_boosting & step <= max.steps) {
            E.hat <- Y - Y.hat
            forest.small <- regression_forest(X,E.hat, sample.fraction = tunable.params["sample.fraction"],
                                        mtry = tunable.params["mtry"], tune.parameters = FALSE,
                                          num.trees = num.trees.tune,
                                           num.threads = num.threads,
                                           min.node.size = tunable.params["min.node.size"],
                                           honesty = honesty,
                                          honesty.fraction = honesty.fraction,
                                           seed = seed, ci.group.size = ci.group.size,
                                          alpha = tunable.params["alpha"],
                                          imbalance.penalty = tunable.params["imbalance.penalty"],
                                        clusters = clusters, samples_per_cluster = samples_per_cluster);
            step.error.approx <- predict(forest.small)$debiased.error
            if((mean(step.error.approx,na.rm=TRUE) < tolerance*mean(error.debiased,na.rm=TRUE))) {
                forest.E <- regression_forest(X,E.hat, sample.fraction = tunable.params["sample.fraction"],
                                            mtry = tunable.params["mtry"], tune.parameters = FALSE,
                                              num.trees = num.trees,
                                               num.threads = num.threads,
                                               min.node.size = tunable.params["min.node.size"],
                                               honesty = honesty,
                                              honesty.fraction = honesty.fraction,
                                               seed = seed, ci.group.size = ci.group.size,
                                              alpha = tunable.params["alpha"],
                                              imbalance.penalty = tunable.params["imbalance.penalty"],
                                            clusters = clusters, samples_per_cluster = samples_per_cluster);
                current.pred <- predict(forest.E)
                Y.hat <- Y.hat + current.pred$predictions
                error.debiased <- current.pred$debiased.error
                boosted.forest[["forests"]][[step]] <- forest.E

            }
            else{
                still_boosting <- FALSE
            }
            step<-step+1
        }
    }

    boosted.forest[["predictions"]] <- Y.hat
    boosted.forest[["debiased.error"]] <- error.debiased

    class(boosted.forest) <- c("boosted_regression_forest")
    boosted.forest
}

#' Predict with a boosted regression forest.
#'
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param Number of boosting iterations to use for prediction. If blank, uses the full number of steps
#'        for the object given
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @return A vector of predictions.
#'
#' @examples \dontrun{
#' # Train a boosted regression forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' r.boosted.forest = boosted_regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' r.pred = predict(r.boosted.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred = predict(r.boosted.forest)
#'
#' }
#'
#' @method predict boosted_regression_forest
#' @export
predict.boosted_regression_forest <- function(object,
                                        newdata=NULL,
                                        num.steps=NULL,
                                      num.threads = NULL
                                      ) {


    # If not on new data, use pre-computed predictions
    if (is.null(newdata)) {
        return(data.frame(predictions=object$predictions
             ,debiased.error=object$debiased.error))
    }
    else {
        forests <- object[["forests"]]
        if (is.null(num.steps)) {
            num.steps <- length(forests)
        }
        num.steps <- min(length(forests),num.steps)
        Y.hat <- predict(forests[[1]],newdata)$predictions
        for (f in 2:num.steps) {
            Y.hat <- Y.hat + predict(forests[[f]],newdata)$predictions
        #    debiased.error <- forest.E$debiased.error
        }
    }

    #not clear how you would get debiased error on new data or why predict returns this
    data.frame(predictions=Y.hat)
}
