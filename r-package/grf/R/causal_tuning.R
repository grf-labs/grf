#' Causal forest tuning
#' 
#' Finds the optimal parameters to be used in training a regression forest. This method
#' currently tunes over min.node.size, mtry, sample.fraction, alpha, and imbalance.penalty.
#' Please see the method 'causal_forest' for a description of the standard causal forest
#' parameters. Note that if fixed values can be supplied for any of the parameters mentioned
#' above, and in that case, that parameter will not be tuned. For example, if this method is
#' called with min.node.size = 10 and alpha = 0.7, then those parameter values will be treated
#' as fixed, and only sample.fraction and imbalance.penalty will be tuned.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param num.fit.reps The number of forests used to fit the tuning model.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction.
#' @param mtry Number of variables tried for each split.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting).
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds 
#'                         to set J1 in the notation of the paper. When using the defaults (honesty = TRUE and 
#'                         honesty.fraction = NULL), half of the data will be used for determining splits
#' @param alpha A tuning parameter that controls the maximum imbalance of a split.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized.
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#'                         determining the imbalance of a split (experimental).
#' @param seed The seed of the C++ random number generator.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#' @param samples_per_cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster. Must be less than the size of the smallest cluster. If set to NULL
#'                            software will set this value to the size of the smallest cluster.#'
#'
#' @return A list consisting of the optimal parameter values ('params') along with their debiased
#'         error ('error').
#'
#' @examples \dontrun{
#' # Find the optimal tuning parameters.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' params = tune_causal_forest(X, Y, W)$params
#'
#' # Use these parameters to train a regression forest.
#' tuned.forest = causal_forest(X, Y, W, num.trees = 1000,
#'     min.node.size = as.numeric(params["min.node.size"]),
#'     sample.fraction = as.numeric(params["sample.fraction"]),
#'     mtry = as.numeric(params["mtry"]),
#'     alpha = as.numeric(params["alpha"]),
#'     imbalance.penalty = as.numeric(params["imbalance.penalty"])
#' }
#'
#' @export
tune_causal_forest <- function(X, Y, W,
                               num.fit.trees = 200,
                               num.fit.reps = 50,
                               num.optimize.reps = 1000,
                               min.node.size = NULL,
                               sample.fraction = 0.5,
                               mtry = NULL,
                               alpha = NULL,
                               imbalance.penalty = NULL,
                               stabilize.splits = TRUE,
                               num.threads = NULL,
                               honesty = TRUE,
                               honesty.fraction = NULL,
                               seed = NULL,
                               clusters = NULL,
                               samples_per_cluster = NULL) {
  validate_X(X)
  if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
  
  num.threads <- validate_num_threads(num.threads)
  seed <- validate_seed(seed)
  clusters <- validate_clusters(clusters, X)
  samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)
  ci.group.size <- 1
  reduced.form.weight <- 0
  honesty.fraction <- validate_honesty_fraction(honesty.fraction, honesty)

  data <- create_data_matrices(X, Y, W)
  outcome.index <- ncol(X) + 1
  treatment.index <- ncol(X) + 2
  instrument.index <- treatment.index
  
  # Separate out the tuning parameters with supplied values, and those that were
  # left as 'NULL'. We will only tune those parameters that the user didn't supply.
  all.params = get_initial_params(min.node.size, sample.fraction, mtry, alpha, imbalance.penalty)
  fixed.params = all.params[!is.na(all.params)]
  tuning.params = all.params[is.na(all.params)]
  
  if (length(tuning.params) == 0) {
    return(list("error"=NA, "params"=c(all.params)))
  }
  
  # Train several mini-forests, and gather their debiased OOB error estimates.
  num.params = length(tuning.params)
  fit.draws = matrix(runif(num.fit.reps * num.params), num.fit.reps, num.params)
  colnames(fit.draws) = names(tuning.params)
  
  debiased.errors = apply(fit.draws, 1, function(draw) {
    params = c(fixed.params, get_params_from_draw(X, draw))
    small.forest <- instrumental_train(data$default, data$sparse,
                                       outcome.index, treatment.index, instrument.index,
                                       as.numeric(params["mtry"]),
                                       num.fit.trees,
                                       num.threads,
                                       as.numeric(params["min.node.size"]),
                                       as.numeric(params["sample.fraction"]),
                                       seed,
                                       honesty,
                                       coerce_honesty_fraction(honesty.fraction),
                                       ci.group.size,
                                       reduced.form.weight,
                                       as.numeric(params["alpha"]),
                                       as.numeric(params["imbalance.penalty"]),
                                       stabilize.splits,
                                       clusters,
                                       samples_per_cluster)
    prediction = instrumental_predict_oob(small.forest, data$default, data$sparse,
                                          num.threads, ci.group.size)
    mean(prediction$debiased.error, na.rm = TRUE)
  })
  
  # Fit the 'dice kriging' model to these error estimates.
  # Note that in the 'km' call, the kriging package prints a large amount of information
  # about the fitting process. Here, capture its console output and discard it.
  variance.guess = rep(var(debiased.errors)/2, nrow(fit.draws))
  env = new.env()
  capture.output(env$kriging.model <-
                   DiceKriging::km(design = data.frame(fit.draws),
                                   response = debiased.errors,
                                   noise.var = variance.guess))
  kriging.model <- env$kriging.model
  
  # To determine the optimal parameter values, predict using the kriging model at a large
  # number of random values, then select those that produced the lowest error.
  optimize.draws = matrix(runif(num.optimize.reps * num.params), num.optimize.reps, num.params)
  colnames(optimize.draws) = names(tuning.params)
  model.surface = predict(kriging.model, newdata=data.frame(optimize.draws), type = "SK")
  
  min.error = min(model.surface$mean)
  optimal.draw = optimize.draws[which.min(model.surface$mean),]
  tuned.params = get_params_from_draw(X, optimal.draw)
  
  list(error = min.error, params = c(fixed.params, tuned.params))
}
