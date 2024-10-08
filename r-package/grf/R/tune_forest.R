#' Tune a forest
#'
#' Finds the optimal parameters to be used in training a forest.
#'
#' @param data The data arguments (output from create_train_matrices) for the forest.
#' @param nrow.X The number of observations.
#' @param ncol.X The number of variables.
#' @param args The remaining call arguments for the forest.
#' @param tune.parameters The vector of parameter names to tune.
#' @param tune.parameters.defaults The grf default values for the vector of parameter names to tune.
#' @param tune.num.trees The number of trees in each 'mini forest' used to fit the tuning model.
#' @param tune.num.reps The number of forests used to fit the tuning model.
#' @param tune.num.draws The number of random parameter values considered when using the model
#'  to select the optimal parameters.
#' @param train The grf forest training function.
#'
#' @return tuning output
#'
#' @keywords internal
tune_forest <- function(data,
                        nrow.X,
                        ncol.X,
                        args,
                        tune.parameters,
                        tune.parameters.defaults,
                        tune.num.trees,
                        tune.num.reps,
                        tune.num.draws,
                        train) {
  fit.parameters <- args[!names(args) %in% tune.parameters]
  fit.parameters[["num.trees"]] <- tune.num.trees
  fit.parameters[["ci.group.size"]] <- 1
  fit.parameters[["compute.oob.predictions"]] <- TRUE
  fit.parameters[["legacy.seed"]] <- get_legacy_seed()

  # 1. Train several mini-forests, and gather their debiased OOB error estimates.
  num.params <- length(tune.parameters)
  unif <- with_seed(runif(tune.num.reps * num.params), seed = args$seed)
  fit.draws <- matrix(unif, tune.num.reps, num.params,
                      dimnames = list(NULL, tune.parameters))

  small.forest.errors <- apply(fit.draws, 1, function(draw) {
    draw.parameters <- get_params_from_draw(nrow.X, ncol.X, draw)
    small.forest <- do.call.rcpp(train, c(data, fit.parameters, draw.parameters))
    error <- small.forest$debiased.error
    mean(error, na.rm = TRUE)
  })
  # Some estimates may be all NaN because of arbitrary inadmissible parameter draws,
  # such as a very low sample.fraction on data with few clusters.
  keep <- !is.na(small.forest.errors)

  if (sum(keep) < 10) {
    warning(paste0(
      "Could not tune forest because nearly all small forest error estimates were NA.\n",
      "Consider increasing tuning arguments tune.num.trees or tune.num.draws."))
    out <- get_tuning_output(params = c(tune.parameters.defaults), status = "failure")
    return(out)
  } else {
    small.forest.errors = small.forest.errors[keep]
    fit.draws = fit.draws[keep, , drop = FALSE]
  }

  if (sd(small.forest.errors) == 0 || sd(small.forest.errors) / mean(small.forest.errors) < 1e-10) {
    warning(paste0(
      "Could not tune forest because small forest errors were nearly constant.\n",
      "Consider increasing argument tune.num.trees."))
    out <- get_tuning_output(params = c(tune.parameters.defaults), status = "failure")
    return(out)
  }

  # 2. Fit the 'dice kriging' model to these error estimates.
  variance.guess <- rep(var(small.forest.errors) / 2, nrow(fit.draws))
  kriging.model <- tryCatch({
    utils::capture.output(
      with_seed(
        model <- DiceKriging::km(
          design = data.frame(fit.draws),
          response = small.forest.errors,
          noise.var = variance.guess
        ),
      seed = args$seed)
    )
    model
  },
  error = function(e) {
    warning(paste0("Dicekriging threw the following error during forest tuning: \n", e))
    return(NULL)
  })

  if (is.null(kriging.model)) {
    warning(paste0("Forest tuning was attempted but failed. Reverting to default parameters."))
    out <- get_tuning_output(params = c(tune.parameters.defaults), status = "failure")
    return(out)
  }

  # 3. To determine the optimal parameter values, predict using the kriging model at a large
  # number of random values, then select those that produced the lowest error.
  unif <- with_seed(runif(tune.num.draws * num.params), seed = args$seed)
  optimize.draws <- matrix(unif, tune.num.draws, num.params,
                           dimnames = list(NULL, tune.parameters))
  model.surface <- predict(kriging.model, newdata = data.frame(optimize.draws), type = "SK")$mean
  tuned.params <- get_params_from_draw(nrow.X, ncol.X, optimize.draws)
  grid <- cbind(error = c(model.surface), tuned.params)
  small.forest.optimal.draw <- which.min(grid[, "error"])

  # To avoid the possibility of selection bias, re-train a moderately-sized forest
  # at the value chosen by the method above
  fit.parameters[["num.trees"]] <- tune.num.trees * 4
  retrained.forest.params <- grid[small.forest.optimal.draw, -1]
  retrained.forest <- do.call.rcpp(train, c(data, fit.parameters, retrained.forest.params))
  retrained.forest.error <- mean(retrained.forest$debiased.error, na.rm = TRUE)

  # 4. Train a forest with default parameters, and check its predicted error.
  # This improves our chances of not doing worse than default
  default.forest <- do.call.rcpp(train, c(data, fit.parameters, tune.parameters.defaults))
  default.forest.error <- mean(default.forest$debiased.error, na.rm = TRUE)

  if (is.na(retrained.forest.error) || default.forest.error < retrained.forest.error) {
    out <- get_tuning_output(
      error = default.forest.error,
      params = tune.parameters.defaults,
      grid = NA,
      status = "default"
    )
  } else {
    out <- get_tuning_output(
      error = retrained.forest.error,
      params = retrained.forest.params,
      grid = grid,
      status = "tuned"
    )
  }

  out
}

get_params_from_draw <- function(nrow.X, ncol.X, draws) {
  if (is.vector(draws)) {
    draws <- rbind(c(draws))
  }
  n <- nrow(draws)
  vapply(colnames(draws), function(param) {
    if (param == "min.node.size") {
      return(floor(2^(draws[, param] * (log(nrow.X) / log(2) - 4))))
    } else if (param == "sample.fraction") {
      return(0.05 + 0.45 * draws[, param])
    } else if (param == "mtry") {
      return(ceiling(min(ncol.X, sqrt(ncol.X) + 20) * draws[, param]))
    } else if (param == "alpha") {
      return(draws[, param] / 4)
    } else if (param == "imbalance.penalty") {
      return(-log(draws[, param]))
    } else if (param == "honesty.fraction") {
      return(0.5 + (0.8 - 0.5) * draws[, param]) # honesty.fraction in U(0.5, 0.8)
    } else if (param == "honesty.prune.leaves") {
      return(ifelse(draws[, param] < 0.5, TRUE, FALSE))
    } else {
      stop("Unrecognized parameter name provided: ", param)
    }
  }, FUN.VALUE = numeric(n))
}

get_tuning_output <- function(status, params, error = NA, grid = NA) {
  out <- list(status = status, params = params, error = error, grid = grid)
  class(out) <- c("tuning_output")
  out
}

with_seed <- function(expr, seed) {
  # Record entry seed
  env <- globalenv()
  oseed <- env$.Random.seed

  # Restore entry seed on exit
  on.exit({
    if (is.null(oseed)) {
      rm(list = ".Random.seed", envir = env, inherits = FALSE)
    } else {
      assign(".Random.seed", value = oseed, envir = env, inherits = FALSE)
    }
  })

  set.seed(seed = seed)
  eval(expr, envir = parent.frame())
}
