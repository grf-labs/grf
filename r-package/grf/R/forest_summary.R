#' Omnibus test for presence of heterogeneity via calibration 
#' 
#' Test calibration of the forest. A slope of 1 means the forest is
#' well calibrated. The significance coefficient of the slope acts as
#' an omnibus test for the presence of heterogeneity: If the slope is
#' not significantly different from 0, then we have not found evidence
#' of heterogeneity.
#'
#' @param forest The trained forest.
#'
#' @return A heteroskedasticity-consistent test of calibration
#'
#' @references Chernozhukov, Victor, Mert Demirer, Esther Duflo, and Ivan Fernandez-Val. "Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized Experiments." arXiv preprint arXiv:1712.04802 (2017).
#' @examples \dontrun{
#' n = 400; p = 5
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.25 + 0.5 * (X[,1] > 0))
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#' forest = causal_forest(X, Y, W)
#' test_calibration(forest)
#' }
#'
#' @export
test_calibration = function(forest) {
    if ("regression_forest" %in% class(forest)) {
        DF = data.frame(target = forest$Y.orig,
                        pred = predict(forest)$predictions)
    } else if ("causal_forest" %in% class(forest)) {
        DF = data.frame(target = forest$Y.orig - forest$Y.hat,
                        pred = (forest$W.orig - forest$W.hat) *
                            predict(forest)$predictions)
        
    }
    
    best.linear.predictor = lm(target ~ pred + 0, data = DF)
    
    # We need to use some extra packages to get
    # heteroskedasticity-consistent inference.
    lmtest::coeftest(best.linear.predictor,
                     vcov = sandwich::vcovHC)
}