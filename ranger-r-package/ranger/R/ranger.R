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
# Institut fuer Medizinische Biometrie und Statistik
# Universitaet zu Luebeck
# Ratzeburger Allee 160
# 23562 Luebeck
# Germany
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
##' The Gini index is used as splitting rule for classification and the estimated response variances for regression.
##' For Survival the log-rank test, a C-index based splitting rule (Schmid et al. 2015) and maximally selected rank statistics (Wright et al. 2016) are available.
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
##' Unordered factor covariates can be handled in 3 different ways by using \code{unordered.factors}: 
##' For 'ignore' all factors are regarded ordered, for 'partition' all possible 2-partitions are considered for splitting and for 'order' the factor levels are ordered by their mean response, as described in Hastie et al. (2009), chapter 9.2.4.
##' The default 'order' is generally recommended, as it computationally fast and can handle an unlimited number of factor levels.
##'
##' For a large number of variables and data frame as input data the formula interface can be slow or impossible to use.
##' Alternatively dependent.variable.name (and status.variable.name for survival) can be used.
##' Consider setting \code{save.memory = TRUE} if you encounter memory problems for very large datasets. 
##' 
##' For GWAS data consider combining \code{ranger} with the \code{GenABEL} package. 
##' See the Examples section below for a demonstration using \code{Plink} data.
##' All SNPs in the \code{GenABEL} object will be used for splitting. 
##' To use only the SNPs without sex or other covariates from the phenotype file, use \code{0} on the right hand side of the formula. 
##' Note that missing values are treated as an extra category while splitting.
##' 
##' See \url{https://github.com/imbs-hl/ranger} for the development version.
##' 
##' To use multithreading on Microsoft Windows platforms, there are currently two options:
##' \itemize{
##'  \item Use R-devel with the new toolchain
##'  \item Install a binary version of ranger, download: \url{https://github.com/imbs-hl/ranger/releases}
##' }
##' 
##' @title Ranger
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param data Training data of class \code{data.frame}, \code{matrix} or \code{gwaa.data} (GenABEL).
##' @param num.trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node. Default is the (rounded down) square root of the number variables. 
##' @param importance Variable importance mode, one of 'none', 'impurity', 'permutation'. The 'impurity' measure is the Gini index for classification and the variance of the responses for regression. For survival, only 'permutation' is available.
##' @param write.forest Save \code{ranger.forest} object, needed for prediction.
##' @param probability Grow a probability forest as in Malley et al. (2012). 
##' @param min.node.size Minimal node size. Default 1 for classification, 5 for regression, 3 for survival, and 10 for probability.
##' @param replace Sample with replacement. 
##' @param sample.fraction Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement. 
##' @param case.weights Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.
##' @param splitrule Splitting rule, survival only. The splitting rule can be chosen of "logrank", "C" and "maxstat" with default "logrank". 
##' @param alpha For "maxstat" splitrule: Significance threshold to allow splitting.
##' @param minprop For "maxstat" splitrule: Lower quantile of covariate distribtuion to be considered for splitting.
##' @param split.select.weights Numeric vector with weights between 0 and 1, representing the probability to select variables for splitting. Alternatively, a list of size num.trees, containing split select weight vectors for each tree can be used.  
##' @param always.split.variables Character vector with variable names to be always tried for splitting.
##' @param unordered.factors Handling of unordered factor covariates, one of 'order', 'partition' and 'ignore' with default 'order'. See below for details. 
##' @param scale.permutation.importance Scale permutation importance by standard error as in (Breiman 2001). Only applicable if permutation variable importance mode selected.
##' @param keep.inbag Save how often observations are in-bag in each tree. 
##' @param holdout Hold-out mode. Hold-out all samples with case weight 0 and use these for variable importance and prediction error.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param save.memory Use memory saving (but slower) splitting mode. No effect for GWAS data.
##' @param verbose Verbose output on or off.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. 
##' @param dependent.variable.name Name of dependent variable, needed if no formula given. For survival forests this is the time variable.
##' @param status.variable.name Name of status variable, only applicable to survival data and needed if no formula given. Use 1 for event and 0 for censoring.
##' @param classification Only needed if data is a matrix. Set to \code{TRUE} to grow a classification forest.
##' @return Object of class \code{ranger} with elements
##'   \item{\code{forest}}{Saved forest (If write.forest set to TRUE). Note that the variable IDs in the \code{split.varIDs} object do not necessarily represent the column number in R.}
##'   \item{\code{predictions}}{Predicted classes/values, based on out of bag samples (classification and regression only).}
##'   \item{\code{forest}}{Saved forest (If write.forest set to TRUE). Note that the variable IDs in the \code{split.varIDs} object do not necessarily represent the column number in R.}
##'   \item{\code{predictions}}{Predicted classes/values, based on out of bag samples (classification and regression only).}
##'   \item{\code{variable.importance}}{Variable importance for each independent variable.}
##'   \item{\code{prediction.error}}{Overall out of bag prediction error. For classification this is the fraction of missclassified samples, for regression the mean squared error and for survival one minus Harrell's c-index.}
##'   \item{\code{r.squared}}{R squared. Also called explained variance or coefficient of determination (regression only).}
##'   \item{\code{confusion.matrix}}{Contingency table for classes and predictions based on out of bag samples (classification only).}
##'   \item{\code{unique.death.times}}{Unique death times (survival only).}
##'   \item{\code{chf}}{Estimated cumulative hazard function for each sample (survival only).}
##'   \item{\code{survival}}{Estimated survival function for each sample (survival only).}
##'   \item{\code{call}}{Function call.}
##'   \item{\code{num.trees}}{Number of trees.}
##'   \item{\code{num.independent.variables}}{Number of independent variables.}
##'   \item{\code{mtry}}{Value of mtry used.}
##'   \item{\code{min.node.size}}{Value of minimal node size used.}
##'   \item{\code{treetype}}{Type of forest/tree. classification, regression or survival.}
##'   \item{\code{importance.mode}}{Importance mode used.}
##'   \item{\code{num.samples}}{Number of samples.}
##'   \item{\code{inbag.counts}}{Number of times the observations are in-bag in the trees.}
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
##' \dontrun{
##' ## Use GenABEL interface to read Plink data into R and grow a classification forest
##' ## The ped and map files are not included
##' library(GenABEL)
##' convert.snp.ped("data.ped", "data.map", "data.raw")
##' dat.gwaa <- load.gwaa.data("data.pheno", "data.raw")
##' phdata(dat.gwaa)$trait <- factor(phdata(dat.gwaa)$trait)
##' ranger(trait ~ ., data = dat.gwaa)
##' }
##'
##' @author Marvin N. Wright
##' @references
##' \itemize{
##'   \item Wright, M. N. & Ziegler, A. (2016). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. Journal of Statistical Software, in press. \url{http://arxiv.org/abs/1508.04409}.
##'   \item Schmid, M., Wright, M. N. & Ziegler, A. (2015). On the Use of Harrell's C for Node Splitting in Random Survival Forests. Technical Report. \url{http://arxiv.org/abs/1507.03092}. 
##'   \item Wright, M. N., Dankowski, T. & Ziegler, A. (2016). Random forests for survival analysis using maximally selected rank statistics. Technical Report. \url{http://arxiv.org/abs/1605.03391}.
##' 
##'   \item Breiman, L. (2001). Random forests. Mach Learn, 45(1), 5-32. 
##'   \item Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. Ann Appl Stat, 841-860. 
##'   \item Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A. (2012). Probability machines: consistent probability estimation using nonparametric learning machines. Methods Inf Med, 51(1), 74.
##'   \item Hastie, T., Tibshirani, R., Friedman, J. (2009). The Elements of Statistical Learning. Springer, New York. 2nd edition.
##'   }
##' @seealso \code{\link{predict.ranger}}
##' @useDynLib ranger
##' @importFrom Rcpp evalCpp
##' @import stats 
##' @import utils
##' @export
ranger <- function(formula = NULL, data = NULL, num.trees = 500, mtry = NULL,
                   importance = "none", write.forest = FALSE, probability = FALSE,
                   min.node.size = NULL, replace = TRUE, 
                   sample.fraction = ifelse(replace, 1, 0.632), 
                   case.weights = NULL, 
                   splitrule = NULL, alpha = 0.5, minprop = 0.1,
                   split.select.weights = NULL, always.split.variables = NULL,
                   unordered.factors = "order",
                   scale.permutation.importance = FALSE,
                   keep.inbag = FALSE, holdout = FALSE,
                   num.threads = NULL, save.memory = FALSE,
                   verbose = TRUE, seed = NULL, 
                   dependent.variable.name = NULL, status.variable.name = NULL, 
                   classification = NULL) {
  
  ## GenABEL GWA data
  if ("gwaa.data" %in% class(data)) {
    snp.names <- data@gtdata@snpnames
    sparse.data <- data@gtdata@gtps@.Data
    data <- data@phdata
    if ("id" %in% names(data)) {
      data$"id" <- NULL
    }
    gwa.mode <- TRUE
    save.memory <- FALSE
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
  
  ## Treetype
  if (is.factor(response)) {
    if (probability) {
      treetype <- 9
    } else {
      treetype <- 1
    }
  } else if (is.numeric(response) & is.vector(response)) {
    if (!is.null(classification) && classification) {
      treetype <- 1
    } else if (probability) {
      treetype <- 9
    } else {
      treetype <- 3
    }
  } else if (class(response) == "Surv" | is.data.frame(response) | is.matrix(response)) {
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
    independent.variable.names <- colnames(data.selected)[colnames(data.selected) != dependent.variable.name &
                                                          colnames(data.selected) != status.variable.name]
  }
  
  ## Recode characters as factors and recode factors if 'order' mode
  if (!is.matrix(data.selected)) {
    character.idx <- sapply(data.selected, is.character)
    
    if (unordered.factors == "order") {
      ## Recode characters and unordered factors
      names.selected <- names(data.selected)
      ordered.idx <- sapply(data.selected, is.ordered)
      factor.idx <- sapply(data.selected, is.factor)
      independent.idx <- names.selected != dependent.variable.name & names.selected != status.variable.name
      recode.idx <- independent.idx & (character.idx | (factor.idx & !ordered.idx))

      ## Recode
      data.selected[recode.idx] <- lapply(data.selected[recode.idx], function(x) {
        if (is.factor(x)) {
          levels <- levels(x)
        } else {
          levels <- unique(x)
        }
        
        ## Order factor levels
        levels.ordered <- names(sort(sapply(levels, function(y) {
          if (treetype == 1 | treetype == 9) {
            mean(as.numeric(response[x == y]))
          } else if (treetype == 3) {
            mean(response[x == y])
          } else if (treetype == 5) {
            mean(response[x == y, 1])
          } 
        })))
        
        ## Return reordered factor
        factor(x, levels = levels.ordered)
      })
      
      ## Save levels
      covariate.levels <- lapply(data.selected[independent.idx], levels)
    } else {
      ## Recode characters only
      data.selected[character.idx] <- lapply(data.selected[character.idx], factor)
    }
  }
  
  ## Input data and variable names
  if (!is.null(formula) & treetype == 5) {
    data.final <- cbind(response[, 1], response[, 2],
                        data.selected[-1])
    colnames(data.final) <- c(dependent.variable.name, status.variable.name,
                              independent.variable.names)
  } else {
    data.final <- data.selected
  }
  if (!is.matrix(data.selected)) {
    ## Create matrix
    data.final <- data.matrix(data.final)
  }
  variable.names <- colnames(data.final)
  
  
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
  
  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  
  ## Keep inbag
  if (!is.logical(keep.inbag)) {
    stop("Error: Invalid value for keep.inbag")
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
  
  ## Sample fraction
  if (!is.numeric(sample.fraction) | sample.fraction <= 0 | sample.fraction > 1) {
    stop("Error: Invalid value for sample.fraction. Please give a value in (0,1].")
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
  
  ## Case weights: NULL for no weights
  if (is.null(case.weights)) {
    case.weights <- c(0,0)
    use.case.weights <- FALSE
    if (holdout) {
      stop("Error: Case weights required to use holdout mode.")
    }
  } else {
    use.case.weights <- TRUE
    
    ## Sample from non-zero weights in holdout mode
    if (holdout) {
      sample.fraction <- sample.fraction * mean(case.weights > 0)
    }
    
    if (!replace && sum(case.weights > 0) < sample.fraction * nrow(data.final)) {
      stop("Error: Fewer non-zero case weights than observations to sample.")
    }
  }
  
  ## Split select weights: NULL for no weights
  if (is.null(split.select.weights)) {
    split.select.weights <- list(c(0,0))
    use.split.select.weights <- FALSE
  } else if (is.numeric(split.select.weights)) {
    if (length(split.select.weights) != length(all.independent.variable.names)) {
      stop("Error: Number of split select weights not equal to number of independent variables.")
    }
    split.select.weights <- list(split.select.weights)
    use.split.select.weights <- TRUE
  } else if (is.list(split.select.weights)) {
    if (length(split.select.weights) != num.trees) {
      stop("Error: Size of split select weights list not equal to number of trees.")
    }
    use.split.select.weights <- TRUE
  } else {
    stop("Error: Invalid split select weights.")
  }
  
  ## Always split variables: NULL for no variables
  if (is.null(always.split.variables)) {
    always.split.variables <- c("0", "0")
    use.always.split.variables <- FALSE
  } else {
    use.always.split.variables <- TRUE
  }
  
  if (use.split.select.weights & use.always.split.variables) {
    stop("Error: Please use only one option of split.select.weights and always.split.variables.")
  }
  
  ## Splitting rule
  if (is.null(splitrule)) {
    splitrule <- "logrank"
    splitrule.num <- 1
  } else if (treetype == 5 & splitrule == "logrank") {
    splitrule.num <- 1
  } else if (treetype == 5 & (splitrule == "auc" | splitrule == "C")) {
    splitrule.num <- 2
  } else if (treetype == 5 & (splitrule == "auc_ignore_ties" | splitrule == "C_ignore_ties")) {
    splitrule.num <- 3
  } else if (treetype == 5 & splitrule == "maxstat") {
    splitrule.num <- 4
  } else {
    stop("Error: Unknown splitrule.")
  }
  
  ## Maxstat splitting
  if (alpha < 0 | alpha > 1) {
    stop("Error: Invalid value for alpha, please give a value between 0 and 1.")
  }
  if (minprop < 0 | minprop > 0.5) {
    stop("Error: Invalid value for minprop, please give a value between 0 and 0.5.")
  }

  ## Unordered factors  
  if (unordered.factors == "partition") {
    names.selected <- names(data.selected)
    ordered.idx <- sapply(data.selected, is.ordered)
    factor.idx <- sapply(data.selected, is.factor)
    independent.idx <- names.selected != dependent.variable.name & names.selected != status.variable.name
    unordered.factor.variables <- names.selected[factor.idx & !ordered.idx & independent.idx]
    
    if (length(unordered.factor.variables) > 0) {
      use.unordered.factor.variables <- TRUE
      ## Check level count
      num.levels <- sapply(data.selected[, factor.idx & !ordered.idx & independent.idx, drop = FALSE], nlevels)
      max.level.count <- 8*.Machine$sizeof.pointer - 1
      if (max(num.levels) > max.level.count) {
        stop(paste("Too many levels in unordered categorical variable ", unordered.factor.variables[which.max(num.levels)], 
                   ". Only ", max.level.count, " levels allowed on this system. Consider using the 'order' option.", sep = ""))
      } 
    } else {
      unordered.factor.variables <- c("0", "0")
      use.unordered.factor.variables <- FALSE
    } 
  } else if (unordered.factors == "ignore" | unordered.factors == "order") {
    ## Ordering for "order" is handled above
    unordered.factor.variables <- c("0", "0")
    use.unordered.factor.variables <- FALSE
  } else {
    stop("Error: Invalid value for unordered.factors, please use 'order', 'partition' or 'ignore'.")
  }
  
  ## Prediction mode always false. Use predict.ranger() method.
  prediction.mode <- FALSE
  predict.all <- FALSE
  
  ## No loaded forest object
  loaded.forest <- list()
  
  ## Clean up
  rm("data.selected")
  
  ## Call Ranger
  result <- rangerCpp(treetype, dependent.variable.name, data.final, variable.names, mtry,
                      num.trees, verbose, seed, num.threads, write.forest, importance.mode,
                      min.node.size, split.select.weights, use.split.select.weights,
                      always.split.variables, use.always.split.variables,
                      status.variable.name, prediction.mode, loaded.forest, sparse.data,
                      replace, probability, unordered.factor.variables, use.unordered.factor.variables, 
                      save.memory, splitrule.num, case.weights, use.case.weights, predict.all, 
                      keep.inbag, sample.fraction, alpha, minprop, holdout)
  
  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }
  
  ## Prepare results
  result$predictions <- drop(do.call(rbind, result$predictions))
  if (importance.mode != 0) {
    names(result$variable.importance) <- all.independent.variable.names
  }
  
  ## Set predictions
  if (treetype == 1 & is.factor(response)) {
    result$predictions <- factor(result$predictions, levels = 1:nlevels(response),
                                 labels = levels(response))
    result$confusion.matrix <- table(unlist(data.final[, dependent.variable.name]), result$predictions, dnn = c("true", "predicted"))
  } else if (treetype == 5) {
    result$chf <- result$predictions
    result$predictions <- NULL
    result$survival <- exp(-result$chf)
    result$splitrule <- splitrule
  } else if (treetype == 9 & !is.matrix(data)) {
    colnames(result$predictions) <- levels(droplevels(response))
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
  result$call <- sys.call()
  result$importance.mode <- importance
  result$num.samples <- nrow(data.final)
  
  ## Write forest object
  if (write.forest) {
    result$forest$levels <- levels(response)
    result$forest$independent.variable.names <- independent.variable.names
    result$forest$treetype <- result$treetype
    class(result$forest) <- "ranger.forest"
    
    ## In 'ordered' mode, save covariate levels
    if (unordered.factors == "order" & !is.matrix(data)) {
      result$forest$covariate.levels <- covariate.levels
    }
  }
  
  class(result) <- "ranger"
  return(result)
}




