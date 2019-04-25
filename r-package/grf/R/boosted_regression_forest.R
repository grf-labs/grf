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
#' @param tolerance If num.steps is NULL, the percentage of previous steps' error that must be estimated
#'                  by cross validation in order to take a new step, default 0.95
#' @param max.steps The maximum number of boosting iterations to try when num.steps NULL
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
                                      tolerance = 0.95,
                                      max.steps = 6,
                                      num.trees.tune = 10) {

  boosted.forest = NULL
  boosted.forest[["forests"]] = list()
  boosted.forest[["debiased.errors"]] = list()
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
  current.pred <- predict(forest.Y)
  #save tuned parameters for use on future boosting iterations
  tunable.params <- forest.Y$tunable.params
  Y.hat <- current.pred$predictions
  error.debiased <- current.pred$debiased.error
  boosted.forest[["forests"]][[1]] <- forest.Y
  boosted.forest[["debiased.errors"]][[1]] <- error.debiased

  still_boosting <- TRUE
  step <- 2

  while(still_boosting) {
    Y.resid <- Y - Y.hat
    #do termination checks
    if(!is.null(num.steps)) {
      still_boosting <- (step <=num.steps)
    } else if(step>max.steps){
      still_boosting <- FALSE
    } else {
      #do cross validation check
      forest.small <- regression_forest(X,Y.resid,
                                        sample.fraction = as.numeric(tunable.params["sample.fraction"]),
                                        mtry = as.numeric(tunable.params["mtry"]), tune.parameters = FALSE,
                                        num.trees = num.trees.tune,
                                        num.threads = num.threads,
                                        min.node.size = as.numeric(tunable.params["min.node.size"]),
                                        honesty = honesty,
                                        honesty.fraction = honesty.fraction,
                                        seed = seed, ci.group.size = ci.group.size,
                                        alpha = as.numeric(tunable.params["alpha"]),
                                        imbalance.penalty = as.numeric(tunable.params["imbalance.penalty"]),
                                        clusters = clusters, samples_per_cluster = samples_per_cluster);
      step.error.approx <- predict(forest.small)$debiased.error
      still_boosting <- mean(step.error.approx,na.rm=TRUE) < tolerance*mean(error.debiased,na.rm=TRUE)
    }

    if(still_boosting) {
      forest.resid <- regression_forest(X,Y.resid, sample.fraction = as.numeric(tunable.params["sample.fraction"]),
                                  mtry = as.numeric(tunable.params["mtry"]), tune.parameters = FALSE,
                                  num.trees = num.trees,
                                  num.threads = num.threads,
                                  min.node.size = as.numeric(tunable.params["min.node.size"]),
                                  honesty = honesty,
                                  honesty.fraction = honesty.fraction,
                                  seed = seed, ci.group.size = ci.group.size,
                                  alpha = as.numeric(tunable.params["alpha"]),
                                  imbalance.penalty = as.numeric(tunable.params["imbalance.penalty"]),
                                  clusters = clusters, samples_per_cluster = samples_per_cluster);

      current.pred <- predict(forest.resid)
      Y.hat <- Y.hat + current.pred$predictions
      error.debiased <- current.pred$debiased.error
      boosted.forest[["forests"]][[step]] <- forest.resid
      boosted.forest[["debiased.errors"]][[step]]<- error.debiased
      step<-step+1
    }
  }
  boosted.forest[["predictions"]] <- Y.hat
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
#' @param num.steps Number of boosting iterations to use for prediction. If blank, uses the full number of steps
#'        for the object given
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
                                              num.steps=NULL) {

  # If not on new data, use pre-computed predictions
  if (is.null(newdata)) {
    return(data.frame(predictions=object$predictions))
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
    }
  }

  data.frame(predictions=Y.hat)
}
