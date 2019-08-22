get_initial_params <- function(min.node.size,
                               sample.fraction,
                               mtry,
                               alpha,
                               imbalance.penalty,
                               honesty,
                               honesty.fraction,
                               prune.empty.leaves) {
  c(
    min.node.size = if (is.null(min.node.size)) NA else validate_min_node_size(min.node.size),
    sample.fraction = if (is.null(sample.fraction)) NA else validate_sample_fraction(sample.fraction),
    mtry = if (is.null(mtry)) NA else validate_mtry(mtry),
    alpha = if (is.null(alpha)) NA else validate_alpha(alpha),
    imbalance.penalty = if (is.null(imbalance.penalty)) NA else validate_imbalance_penalty(imbalance.penalty),
    honesty.fraction = if (is.null(honesty.fraction) && honesty) NA else
     validate_honesty_fraction(honesty.fraction, honesty),
    prune.empty.leaves = if (is.null(prune.empty.leaves) && honesty) NA else validate_prune_empty_leaves(prune.empty.leaves)
  )
}

get_params_from_draw <- function(X, draws) {
  if (is.vector(draws)) {
    draws <- rbind(c(draws))
  }
  n <- nrow(draws)
  vapply(colnames(draws), function(param) {
    if (param == "min.node.size") {
      return(floor(2^(draws[, param] * (log(nrow(X)) / log(2) - 4))))
    } else if (param == "sample.fraction") {
      return(0.05 + 0.45 * draws[, param])
    } else if (param == "mtry") {
      return(ceiling(min(ncol(X), sqrt(ncol(X)) + 20) * draws[, param]))
    } else if (param == "alpha") {
      return(draws[, param] / 4)
    } else if (param == "imbalance.penalty") {
      return(-log(draws[, param]))
    } else if (param == "honesty.fraction") {
      return(0.5 + (0.8 - 0.5) * draws[, param]) # honesty.fraction in U(0.5, 0.8)
    } else if (param == "prune.empty.leaves") {
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
