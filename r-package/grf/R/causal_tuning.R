#' Causal forest tuning
#' 
#' Finds the optimal parameters to be used in training a regression forest. This method
#' currently tunes over min.node.size, sample.fraction, alpha, and imbalance.penalty.
#' Please see the method 'regression_forest' for a description of the standard forest
#' parameters. Note that if fixed values can be supplied for any of the parameters mentioned
#' above, and in that case, that parameter will not be tuned. For example, if this method is
#' called with min.node.size = 10 and alpha = 0.7, then those parameter values will be treated
#' as fixed, and only sample.fraction and imbalance.penalty will be tuned.
#'
#' @param num.fit.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param num.fit.reps The number of forests used to fit the tuning model.
#' @param num.optimize.reps The number of random parameter values considered when using the model
#'                          to select the optimal parameters.
#'
#' @return A list consisting of the optimal parameter values ('params') along with their debiased
#'         error ('error').
#'
#' @examples \dontrun{
#' # Find the optimal tuning parameters.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' params = tune_regression_forest(X, Y)$params
#'
#' # Use these parameters to train a regression forest.
#' tuned.forest = causal_forest(X, Y, num.trees = 1000,
#'     min.node.size = params["min.node.size"],
#'     sample.fraction = params["sample.fraction"],
#'     alpha = params["alpha"],
#'     alpha = params["alpha"],
#'     imbalance.penalty = params["imbalance.penalty"])
#' }
#'
#' @export
tune_causal_forest <- function(X, Y, W,
                               num.fit.trees = 40,
                               num.fit.reps = 100,
                               num.optimize.reps = 1000,
                               min.node.size = NULL,
                               sample.fraction = NULL,
                               mtry = NULL,
                               alpha = NULL,
                               imbalance.penalty = NULL,
                               stabilize.splits = TRUE,
                               num.threads = NULL,
                               honesty = TRUE,
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

  #print(debiased.errors)
  
  # Fit the 'dice kriging' model to these error estimates.
  # Note that in the 'km' call, the kriging package prints a large amount of information
  # about the fitting process. Here, capture its console output and discard it.
  variance.guess = rep(var(debiased.errors)/2, nrow(fit.draws))
  env = environment()
  capture.output(env$kriging.model <-
                   DiceKriging::km(design = data.frame(fit.draws),
                                   response = debiased.errors,
                                   noise.var = variance.guess))
  
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
