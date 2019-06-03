#' Omnibus evaluation of the quality of the random forest estimates via calibration.
#'
#' Test calibration of the forest. Computes the best linear fit of the target
#' estimand using the forest prediction (on held-out data) as well as the mean
#' forest prediction as the sole two regressors. A coefficient of 1 for
#' `mean.forest.prediction` suggests that the mean forest prediction is correct,
#' whereas a coefficient of 1 for `differential.forest.prediction` additionally suggests
#' that the forest has captured heterogeneity in the underlying signal.
#' The p-value of the `differential.forest.prediction` coefficient
#' also acts as an omnibus test for the presence of heterogeneity: If the coefficient
#' is significantly greater than 0, then we can reject the null of
#' no heterogeneity.
#'
#' @param forest The trained forest.
#' @return A heteroskedasticity-consistent test of calibration.
#' @references Chernozhukov, Victor, Mert Demirer, Esther Duflo, and Ivan Fernandez-Val.
#'             "Generic Machine Learning Inference on Heterogenous Treatment Effects in
#'             Randomized Experiments." arXiv preprint arXiv:1712.04802 (2017).
#'
#' @examples
#' \dontrun{
#' n <- 800
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.25 + 0.5 * (X[, 1] > 0))
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' forest <- causal_forest(X, Y, W)
#' test_calibration(forest)
#' }
#'
#' @export
test_calibration <- function(forest) {
  observation.weight <- observation_weights(forest)
  clusters <- if (length(forest$clusters) > 0) {
    forest$clusters
  } else {
    1:length(observation.weight)
  }
  if ("regression_forest" %in% class(forest)) {
    preds <- predict(forest)$predictions
    mean.pred <- weighted.mean(preds, observation.weight)
    DF <- data.frame(
      target = unname(forest$Y.orig),
      mean.forest.prediction = mean.pred,
      differential.forest.prediction = preds - mean.pred
    )
  } else if ("causal_forest" %in% class(forest)) {
    preds <- predict(forest)$predictions
    mean.pred <- weighted.mean(preds, observation.weight)
    DF <- data.frame(
      target = unname(forest$Y.orig - forest$Y.hat),
      mean.forest.prediction = unname(forest$W.orig - forest$W.hat) * mean.pred,
      differential.forest.prediction = unname(forest$W.orig - forest$W.hat) *
        (preds - mean.pred)
    )
  } else {
    stop("Calibration check not supported for this type of forest.")
  }

  best.linear.predictor <-
    lm(target ~ mean.forest.prediction + differential.forest.prediction + 0,
      weights = observation.weight,
      data = DF
    )
  blp.summary <- lmtest::coeftest(best.linear.predictor,
    vcov = sandwich::vcovCL,
    type = "HC3",
    cluster = clusters
  )
  attr(blp.summary, "method") <-
    paste("Best linear fit using forest predictions (on held-out data)",
      "as well as the mean forest prediction as regressors, along",
      "with one-sided heteroskedasticity-robust (HC3) SEs",
      sep = "\n"
    )
  # convert to one-sided p-values
  dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
  blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[, 4] / 2, blp.summary[, 4] / 2)
  blp.summary
}
