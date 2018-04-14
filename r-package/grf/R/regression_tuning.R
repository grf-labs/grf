#' Regression forest tuning
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
#' tuned.forest = regression_forest(X, Y, num.trees = 1000,
#'     min.node.size = params["min.node.size"],
#'     sample.fraction = params["sample.fraction"],
#'     alpha = params["alpha"],
#'     alpha = params["alpha"],
#'     imbalance.penalty = params["imbalance.penalty"])
#' }
#'
#' @export
tune_regression_forest <- function(X, Y,
                                   num.fit.trees = 10,
                                   num.fit.reps = 100,
                                   num.optimize.reps = 1000,
                                   min.node.size = NULL,
                                   sample.fraction = NULL,
                                   mtry = NULL,
                                   alpha = NULL,
                                   imbalance.penalty = NULL,
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
  
  data <- create_data_matrices(X, Y)
  outcome.index <- ncol(X) + 1
  
  # Separate out the tuning parameters with supplied values, and those that were
  # left as 'NULL'. We will only tune those parameters that the user didn't supply.
  all.params = c(min.node.size = validate_min_node_size(min.node.size, null_to_na = TRUE),
                 sample.fraction = validate_sample_fraction(sample.fraction, null_to_na = TRUE),
                 mtry = validate_mtry(mtry, X, null_to_na = TRUE),
                 alpha = validate_alpha(alpha, null_to_na = TRUE),
                 imbalance.penalty = validate_imbalance_penalty(imbalance.penalty, null_to_na = TRUE))
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
    params = c(fixed.params, get_params(X, draw))
    small.forest = regression_forest(X, Y, tune.parameters = FALSE,
                                     num.threads = num.threads, honesty = honesty, seed = seed,
                                     clusters = clusters, samples_per_cluster = samples_per_cluster,
                                     num.trees = num.fit.trees,
                                     min.node.size = as.numeric(params["min.node.size"]), 
                                     sample.fraction = as.numeric(params["sample.fraction"]),
                                     mtry = as.numeric(params["mtry"]),
                                     alpha = as.numeric(params["alpha"]),
                                     imbalance.penalty = as.numeric(params["imbalance.penalty"])
    )
    prediction = predict(small.forest)
    mean(prediction$debiased.error, na.rm = TRUE)
  })
  
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
  fitted.params = get_params(X, optimal.draw)
  
  list("error"=min.error, "params"=c(fixed.params, fitted.params))
}


get_params <- function(X, draw) {
  result = c()
  for (param in names(draw)) {
    if (param == "min.node.size") {
      value = floor(2^(draw[param] * (log(nrow(X)) / log(2) - 4)))
    } else if (param == "sample.fraction") {
      value = 0.05 + 0.45 * draw[param]
    } else if (param == "mtry") {
      value = ceiling(ncol(X) * draw[param])
    } else if (param == "alpha") {
      value = draw[param]/4
    } else if (param == "imbalance.penalty") {
      value = -log(draw[param])
    } else {
      stop("Unrecognized parameter name provided: ", param)
    }
    result = c(result, value)
  }
  result
}