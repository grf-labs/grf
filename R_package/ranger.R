# -------------------------------------------------------------------------------
#   This file is part of Ranger.
#
# Ranger is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ranger is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ranger. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Marvin N. Wright
# Institut für Medizinische Biometrie und Statistik
# Universität zu Lübeck
# Ratzeburger Allee 160
# 23562 Lübeck
#
# http://www.imbs-luebeck.de
# wright@imbs.uni-luebeck.de
# -------------------------------------------------------------------------------

##' Ranger is a fast implementation of Random Forest (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data.
##' Classification, regression, and survival forests are supported.
##' Classification and regression forests are implemented as in the original Random Forest (Breiman 2001), survival forests as in Random Survival Forests (Ishwaran et al. 2008).
##'
##' The tree type is determined by the type of the dependent variable.
##' For factors classification trees are grown, for numeric values regression trees and for survival objects survival trees.
##' The Gini index is used as splitting rule for classification, the estimated response variances for regression and the log-rank test for survival.
##' For Survival the log-rank test or an AUC-based splitting rule are available.
##'
##' With the \code{probability} option and factor dependent variable a probability forest is grown.
##' Here, the estimated response variances are used for splitting, as in regression forests.
##' Predictions are class probabilities for each sample.
##' For details see Malley et al. (2012).
##'
##' Note that for classification and regression nodes with size smaller than min.node.size can occur, like in original Random Forest.
##' For survival all nodes contain at least min.node.size samples. 
##' Variables selected with \code{always.split.variables} are tried additionaly to the mtry variables randomly selected.
##' In \code{split.select.weights} variables weighted with 0 are never selected and variables with 1 are always selected. 
##' Weights do not need to sum up to 1, they will be normalized later. 
##' The usage of \code{split.select.weights} can increase the computation times for large forests.
##'
##' For a large number of variables and data frame as input data the formula interface can be slow.
##' Alternatively dependent.variable.name (and status.variable.name for survival) can be used.
##' By setting \code{memory} to 'float' or 'char' the memory usage can be reduced with a possible loss of precision. 
##' Use char only if all endpoints and covariates are integer values between -128 and 127 or factors with few levels. 
##' 
##' Multithreading is currently not supported for Microsoft Windows platforms. 
##'
##' @title Ranger
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param data Training data of class \code{data.frame} or \code{gwaa.data} (GenABEL).
##' @param num.trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node.
##' @param importance Variable importance mode, one of 'none', 'impurity', 'permutation'. The 'impurity' measure is the Gini index for classification and the variance of the responses for regression.
##' @param write.forest Save \code{ranger.forest} object, needed for prediction.
##' @param probability Grow a probability forest. This is a classification forest which returns class probabilities instead of classifications.
##' @param min.node.size Minimal node size. Default 1 for classification, 5 for regression, 3 for survival, and 10 for probability.
##' @param replace Sample with replacement. Default TRUE.
##' @param splitrule Splitting rule, survival only. The splitting rule can be chosen of "logrank", "auc" and "auc_ignore_ties" with default "logrank". 
##' @param split.select.weights Numeric vector with weights between 0 and 1, representing the probability to select variables for splitting.  
##' @param always.split.variables Character vector with variable names to be always tried for splitting.
##' @param scale.permutation.importance Scale permutation importance by standard error as in (Breiman 2001). Only applicable if permutation variable importance mode selected.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param verbose Verbose output on or off.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. 
##' @param memory Memory mode, one of 'double', 'float', 'char'. Default 'double'.
##' @param dependent.variable.name Name of dependent variable, needed if no formula given. For survival forests this is the time variable.
##' @param status.variable.name Name of status variable, only applicable to survival data and needed if no formula given. Use 1 for event and 0 for censoring.
##' @return Object of class \code{ranger} with elements
##'   \tabular{ll}{
##'       \code{forest} \tab Saved forest (If write.forest set to TRUE). \cr
##'       \code{predictions}    \tab Predicted classes/values, based on out of bag samples (classification and regression only). \cr
##'       \code{variable.importance}     \tab Variable importance for each independent variable. \cr
##'       \code{prediction.error}   \tab Overall out of bag prediction error. For classification this is the fraction of missclassified samples, for regression the mean squared error and for survival one minus Harrell's c-index. \cr
##'       \code{r.squared}   \tab R squared. Also called explained variance or coefficient of determination (regression only). \cr
##'       \code{classification.table} \tab Contingency table for classes and predictions based on out of bag samples (classification only). \cr
##'       \code{unique.death.times} \tab Unique death times (survival only). \cr
##'       \code{chf} \tab Estimated cumulative hazard function for each sample (survival only). \cr
##'       \code{survival} \tab Estimated survival function for each sample (survival only). \cr
##'       \code{call}    \tab Function call. \cr
##'       \code{num.trees}   \tab Number of trees. \cr
##'       \code{num.independent.variables} \tab Number of independent variables. \cr
##'       \code{mtry}    \tab Value of mtry used. \cr
##'       \code{min.node.size}   \tab Value of minimal node size used. \cr
##'       \code{treetype}    \tab Type of forest/tree. classification, regression or survival. \cr
##'       \code{memory.mode}   \tab Memory mode used. \cr
##'       \code{importance.mode}     \tab Importance mode used. \cr
##'       \code{num.samples}     \tab Number of samples.
##'   }
##' @examples
##' require(ranger)
##'
##' ## Classification forest with default settings
##' ranger(Species ~ ., data = iris)
##'
##' ## Prediction
##' train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
##' iris.train <- iris[train.idx, ]
##' iris.test <- iris[-train.idx, ]
##' rg.iris <- ranger(Species ~ ., data = iris.train, write.forest = TRUE)
##' pred.iris <- predict(rg.iris, dat = iris.test)
##' table(iris.test$Species, pred.iris$predictions)
##'
##' ## Variable importance
##' rg.iris <- ranger(Species ~ ., data = iris, importance = "impurity")
##' rg.iris$variable.importance
##'
##' ## Survival forest
##' require(survival)
##' rg.veteran <- ranger(Surv(time, status) ~ ., data = veteran)
##' plot(rg.veteran$unique.death.times, rg.veteran$survival[1,])
##'
##' ## Alternative interface
##' ranger(dependent.variable.name = "Species", data = iris)
##'
##' @author Marvin N. Wright
##' @references
##'   Breiman, L. (2001). Random forests. Mach Learn, 45(1), 5-32. \cr
##'   Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. Ann Appl Stat, 841-860. \cr
##'   Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A. (2012). Probability machines: consistent probability estimation using nonparametric learning machines. Methods Inf Med, 51(1), 74.
##' @seealso \code{\link{predict.ranger}}
##' @export
ranger <- function(formula = NULL, data = NULL, num.trees = 500, mtry = NULL,
                    importance = "none", write.forest = FALSE, probability = FALSE,
                    min.node.size = NULL, replace = TRUE, splitrule = NULL,
                    split.select.weights = NULL, always.split.variables = NULL,
                    scale.permutation.importance = FALSE,
                    num.threads = NULL,
                    verbose = TRUE, seed = NULL, memory = "double",
                    dependent.variable.name = NULL, status.variable.name = NULL) {
  
  ## GenABEL GWA data
  if (class(data) == "gwaa.data") {
    snp.names <- data@gtdata@snpnames
    sparse.data <- data@gtdata@gtps@.Data
    data <- data@phdata
    if ("id" %in% names(data)) {
      data$"id" <- NULL
    }
    gwa.mode <- TRUE
  } else {
    sparse.data <- as.matrix(0)
    gwa.mode <- FALSE
  }
  
  ## Formula interface. Use whole data frame is no formula provided and depvarname given
  if (is.null(formula)) {
    if (is.null(dependent.variable.name)) {
      stop("Error: Please give formula or dependent variable name.")
    }
    if (is.null(status.variable.name)) {
      status.variable.name <- "none"
      response <- data[, dependent.variable.name]
    } else {
      response <- data[, c(dependent.variable.name, status.variable.name)]
    }
    data.selected <- data
  } else {
    formula <- formula(formula)
    if (class(formula) != "formula") {
      stop("Error: Invalid formula.")
    }
    data.selected <- model.frame(formula, data, na.action = na.fail)
    response <- data.selected[[1]]
  }
  
  ## Probability estimation
  if (probability & !is.factor(response)) {
    stop("Error: Probability estimation is only applicable to categorical (factor) dependent variables.")
  }
  
  ## Treetype
  if (is.factor(response)) {
    if (probability) {
      treetype <- 9
    } else {
      treetype <- 1
    }
  } else if (is.numeric(response) & is.vector(response)) {
    treetype <- 3
  } else if (class(response) == "Surv" | class(response) == "data.frame") {
    treetype <- 5
  } else {
    stop("Error: Unsupported type of dependent variable.")
  }
  
  ## Dependent and status variable name. For non-survival dummy status variable name.
  if (!is.null(formula)) {
    if (treetype == 5) {
      dependent.variable.name <- dimnames(response)[[2]][1]
      status.variable.name <- dimnames(response)[[2]][2]
    } else {
      dependent.variable.name <- names(data.selected)[1]
      status.variable.name <- "none"
    }
    independent.variable.names <- names(data.selected)[-1]
  } else {
    independent.variable.names <- names(data.selected)[names(data.selected) != dependent.variable.name &
                                                         names(data.selected) != status.variable.name]
  }
  
  ## Input data and variable names
  if (!is.null(formula)) {
    if (treetype == 5) {
      data.final <- data.matrix(cbind(response[, 1], response[, 2],
                                      data.selected[-1]))
      variable.names <- c(dependent.variable.name, status.variable.name,
                          independent.variable.names)
    } else {
      data.final <- data.matrix(data.selected)
      variable.names <- names(data.selected)
    }
  } else {
    data.final <- data.matrix(data.selected)
    variable.names <- names(data.selected)
  }
  
  ## If gwa mode, add snp variable names
  if (gwa.mode) {
    variable.names <- c(variable.names, snp.names)
    all.independent.variable.names <- c(independent.variable.names, snp.names)
  } else {
    all.independent.variable.names <- independent.variable.names
  }
  
  ## Number of trees
  if (!is.numeric(num.trees) | num.trees < 1) {
    stop("Error: Invalid value for num.trees.")
  }
  
  ## mtry
  if (is.null(mtry)) {
    mtry <- 0
  } else if (!is.numeric(mtry) | mtry < 0) {
    stop("Error: Invalid value for mtry")
  }
  
  ## Memory mode
  if (is.null(memory) | memory == "double") {
    memory.mode <- 0
  } else if (memory == "float") {
    memory.mode <- 1
  } else if (memory == "char") {
    memory.mode <- 2
  } else {
    stop("Error: Unknown memory mode.")
  }
  
  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  
  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) | num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }
  
  ## Minumum node size
  if (is.null(min.node.size)) {
    min.node.size <- 0
  } else if (!is.numeric(min.node.size) | min.node.size < 0) {
    stop("Error: Invalid value for min.node.size")
  }
  
  ## Importance mode
  if (is.null(importance) | importance == "none") {
    importance.mode <- 0
  } else if (importance == "impurity") {
    importance.mode <- 1
    if (treetype == 5) {
      stop("Node impurity variable importance not supported for survival forests.")
    }
  } else if (importance == "permutation") {
    if (scale.permutation.importance) {
      importance.mode <- 2
    } else {
      importance.mode <- 3
    }
  } else {
    stop("Error: Unknown importance mode.")
  }
  
  ## Split select weights: NULL for no weights
  if (is.null(split.select.weights)) {
    split.select.weights <- c(0,0)
    use.split.select.weights <- FALSE
  } else {
    use.split.select.weights <- TRUE
  }
  
  ## Always split variables: NULL for no variables
  if (is.null(always.split.variables)) {
    always.split.variables <- c("0", "0")
    use.always.split.variables <- FALSE
  } else {
    use.always.split.variables <- TRUE
  }
  
  if (use.split.select.weights & use.always.split.variables) {
    stop("Error: Please use only one option of use.split.select.weights and use.always.split.variables.")
  }
  
  ## Splitting rule
  if (is.null(splitrule)) {
    splitrule <- 1
  } else if (treetype == 5 & splitrule == "logrank") {
    splitrule <- 1
  } else if (treetype == 5 & splitrule == "auc") {
    splitrule <- 2
  } else if (treetype == 5 & splitrule == "auc_ignore_ties") {
    splitrule <- 3
  } else {
    stop("Error: Unknown splitrule.")
  }

  ## Prediction mode always false. Use predict.ranger() method.
  prediction.mode <- FALSE
  
  ## No loaded forest object
  loaded.forest <- list()
  
  ## Call Ranger
  result <- rangerCpp(treetype, dependent.variable.name, memory.mode, data.final, variable.names, mtry,
              num.trees, verbose, seed, num.threads, write.forest, importance.mode,
              min.node.size, split.select.weights, use.split.select.weights,
              always.split.variables, use.always.split.variables,
              status.variable.name, prediction.mode, loaded.forest, sparse.data,
              replace, probability, splitrule)

  if (length(result) == 0) {
    stop("Internal error.")
  }
  
  ## Prepare results
  result$predictions <- drop(do.call(rbind, result$predictions))
  if (importance.mode != 0) {
    names(result$variable.importance) <- all.independent.variable.names
  }
  
  ## Set predictions
  if (treetype == 1) {
    result$predictions <- factor(result$predictions, levels = 1:nlevels(response),
                                 labels = levels(response))
    result$classification.table <- table(result$predictions, unlist(data[, dependent.variable.name]), dnn = c("predicted", "true"))
  } else if (treetype == 5) {
    result$chf <- result$predictions
    result$predictions <- NULL
    result$survival <- exp(-result$chf)
  } else if (treetype == 9) {
    colnames(result$predictions) <- levels(response)
  }
  
  ## Set treetype
  if (treetype == 1) {
    result$treetype <- "Classification"
  } else if (treetype == 3) {
    result$treetype <- "Regression"
  } else if (treetype == 5) {
    result$treetype <- "Survival"
  } else if (treetype == 9) {
    result$treetype <- "Probability estimation"
  }
  if (treetype == 3) {
    result$r.squared <- 1 - result$prediction.error / var(response)
  }
  result$call <- match.call()
  result$memory.mode <- memory
  result$importance.mode <- importance
  result$num.samples <- nrow(data.final)
  
  ## Write forest object
  if (write.forest) {
    result$forest$levels <- levels(response)
    result$forest$independent.variable.names <- independent.variable.names
    result$forest$treetype <- result$treetype
    class(result$forest) <- "ranger.forest"
  }
  
  class(result) <- "ranger"
  return(result)
}




