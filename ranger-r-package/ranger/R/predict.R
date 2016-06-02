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

##' Prediction with new data and a saved forest from Ranger.
##'
##' @title Ranger prediction
##' @param object Ranger \code{ranger.forest} object.
##' @param data New test data of class \code{data.frame} or \code{gwaa.data} (GenABEL).
##' @param predict.all Return a matrix with individual predictions for each tree instead of aggregated predictions for all trees (classification and regression only).
##' @param seed Random seed used in Ranger.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param verbose Verbose output on or off.
##' @param ... further arguments passed to or from other methods.
##' @return Object of class \code{ranger.prediction} with elements
##'   \tabular{ll}{
##'       \code{predictions}    \tab Predicted classes/values (only for classification and regression)  \cr
##'       \code{unique.death.times} \tab Unique death times (only for survival). \cr
##'       \code{chf} \tab Estimated cumulative hazard function for each sample (only for survival). \cr
##'       \code{survival} \tab Estimated survival function for each sample (only for survival). \cr
##'       \code{num.trees}   \tab Number of trees. \cr
##'       \code{num.independent.variables} \tab Number of independent variables. \cr
##'       \code{treetype}    \tab Type of forest/tree. Classification, regression or survival. \cr
##'       \code{num.samples}     \tab Number of samples.
##'   }
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @export
predict.ranger.forest <- function(object, data, predict.all = FALSE,
                                  seed = NULL, num.threads = NULL,
                                  verbose = TRUE, ...) {
  
  ## GenABEL GWA data
  if ("gwaa.data" %in% class(data)) {
    snp.names <- snp.names(data)
    sparse.data <- data@gtdata@gtps@.Data
    data <- data@phdata[, -1]
    gwa.mode <- TRUE
    variable.names <- c(names(data), snp.names)
  } else {
    sparse.data <- as.matrix(0)
    gwa.mode <- FALSE
    variable.names <- colnames(data)
  }

  ## Check forest argument
  if (class(object) != "ranger.forest") {
    stop("Error: Invalid class of input object.")
  } else {
    forest <- object
  }
  if (is.null(forest$dependent.varID) | is.null(forest$num.trees) |
        is.null(forest$child.nodeIDs)  | is.null(forest$split.varIDs) |
        is.null(forest$split.values) | is.null(forest$independent.variable.names) |
        is.null(forest$treetype)) {
    stop("Error: Invalid forest object.")
  }
  if (forest$treetype == "Survival" & (is.null(forest$status.varID)  |
                                         is.null(forest$chf) | is.null(forest$unique.death.times))) {
    stop("Error: Invalid forest object.")
  }
  
  ## Create final data
  if (forest$treetype == "Survival") {
    if (forest$dependent.varID > 0 & forest$status.varID > 1) {
      ## If alternative interface used, don't subset data
      data.used <- data
    } else {
      ## If formula interface used, subset data
      data.selected <- subset(data, select = forest$independent.variable.names)
      
      ## Arange data as in original data
      data.used <- cbind(0, 0, data.selected)
      variable.names <- c("time", "status", forest$independent.variable.names)
    }
      
  ## Index of no-recode variables
  idx.norecode <- c(-(forest$dependent.varID+1), -(forest$status.varID+1))

  } else {
    ## No survival
    if (ncol(data) == length(forest$independent.variable.names)+1 & forest$dependent.varID > 0) {
      ## If alternative interface used, don't subset data
      data.used <- data
    } else {
      ## If formula interface used, subset data
      data.selected <- subset(data, select = forest$independent.variable.names)
      
      ## Arange data as in original data
      if (forest$dependent.varID == 0) {
        data.used <- cbind(0, data.selected)
        variable.names <- c("dependent", forest$independent.variable.names)
      } else if (forest$dependent.varID >= ncol(data)) {
        data.used <- cbind(data.selected, 0)
        variable.names <- c(forest$independent.variable.names, "dependent")
      } else {
        data.used <- cbind(data.selected[, 1:forest$dependent.varID], 
                           0, 
                           data.selected[, (forest$dependent.varID+2):ncol(data.selected)])
        variable.names <- c(forest$independent.variable.names[1:forest$dependent.varID], 
                            "dependent", 
                            forest$independent.variable.names[(forest$dependent.varID+1):length(forest$independent.variable.names)])
      }
    }
    
    ## Index of no-recode variables
    idx.norecode <- -(forest$dependent.varID+1)
  }
  
  ## Recode characters
  if (!is.matrix(data.used)) {
    char.columns <- sapply(data.used, is.character)
    data.used[char.columns] <- lapply(data.used[char.columns], factor)
  }

  ## Recode factors if forest grown 'order' mode
  if (!is.null(forest$covariate.levels) && !all(sapply(forest$covariate.levels, is.null))) {
    data.used[, idx.norecode] <- mapply(function(x, y) {
      if(is.null(y)) {
        x
      } else {
        factor(x, levels = y)
      }
    }, data.used[, idx.norecode], forest$covariate.levels, SIMPLIFY = FALSE)
  }
  
  ## Convert to data matrix
  data.final <- data.matrix(data.used)

  ## If gwa mode, add snp variable names
  if (gwa.mode) {
    variable.names <- c(variable.names, snp.names)
  }

  if (any(is.na(data.final))) {
    stop("Missing values in data.")
  }
  
  if (sum(!(forest$independent.variable.names %in% variable.names)) > 0) {
    stop("Error: One or more independent variables not found in data.")
  }

  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) | num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }

  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }

  if (forest$treetype == "Classification") {
    treetype <- 1
  } else if (forest$treetype == "Regression") {
    treetype <- 3
  } else if (forest$treetype == "Survival") {
    treetype <- 5
  } else if (forest$treetype == "Probability estimation") {
    treetype <- 9
  } else {
    stop("Error: Unknown tree type.")
  }

  ## Defaults for variables not needed
  dependent.variable.name <- "none"
  mtry <- 0
  importance <- 0
  min.node.size <- 0
  split.select.weights <- list(c(0, 0))
  use.split.select.weights <- FALSE
  always.split.variables <- c("0", "0")
  use.always.split.variables <- FALSE
  status.variable.name <- "status"
  prediction.mode <- TRUE
  write.forest <- FALSE
  replace <- TRUE
  probability <- FALSE
  unordered.factor.variables <- c("0", "0")
  use.unordered.factor.variables <- FALSE
  save.memory <- FALSE
  splitrule <- 1
  alpha <- 0
  minprop <- 0
  case.weights <- c(0, 0)
  use.case.weights <- FALSE
  keep.inbag <- FALSE
  sample.fraction <- 1
  holdout <- FALSE
  
  ## Call Ranger
  result <- rangerCpp(treetype, dependent.variable.name, data.final, variable.names, mtry,
                      forest$num.trees, verbose, seed, num.threads, write.forest, importance,
                      min.node.size, split.select.weights, use.split.select.weights,
                      always.split.variables, use.always.split.variables,
                      status.variable.name, prediction.mode, forest, sparse.data, replace, probability,
                      unordered.factor.variables, use.unordered.factor.variables, save.memory, splitrule, 
                      case.weights, use.case.weights, predict.all, keep.inbag, sample.fraction, 
                      alpha, minprop, holdout)

  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }
  
  ## Prepare results
  result$predictions <- drop(do.call(rbind, result$predictions))
  result$num.samples <- nrow(data.final)
  result$treetype <- forest$treetype
  
  if (forest$treetype == "Classification" & !is.null(forest$levels)) {
    if (!predict.all) {
      result$predictions <- factor(result$predictions, levels = 1:length(forest$levels),
                                   labels = forest$levels)
    }
  } else if (forest$treetype == "Survival") {
    result$unique.death.times <- forest$unique.death.times
    result$chf <- result$predictions
    result$predictions <- NULL
    result$survival <- exp(-result$chf)
  } else if (forest$treetype == "Probability estimation" & !is.null(forest$levels)) {
    if (is.matrix(result$predictions)) {
      colnames(result$predictions) <- forest$levels
    } else {
      names(result$predictions) <- forest$levels
    }
    
  }

  class(result) <- "ranger.prediction"
  return(result)
}

##' Prediction with new data and a saved forest from Ranger.
##' 
##' For classification and predict.all = TRUE, a matrix of factor levels is returned. 
##' To retrieve the corresponding factor levels, use rf$forest$levels, if rf is the ranger object.
##'
##' @title Ranger prediction
##' @param object Ranger \code{ranger} object.
##' @param data New test data of class \code{data.frame} or \code{gwaa.data} (GenABEL).
##' @param predict.all Return a matrix with individual predictions for each tree instead of aggregated predictions for all trees (classification and regression only).
##' @param seed Random seed used in Ranger.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param verbose Verbose output on or off.
##' @param ... further arguments passed to or from other methods.
##' @return Object of class \code{ranger.prediction} with elements
##'   \tabular{ll}{
##'       \code{predictions}    \tab Predicted classes/values (only for classification and regression)  \cr
##'       \code{unique.death.times} \tab Unique death times (only for survival). \cr
##'       \code{chf} \tab Estimated cumulative hazard function for each sample (only for survival). \cr
##'       \code{survival} \tab Estimated survival function for each sample (only for survival). \cr
##'       \code{num.trees}   \tab Number of trees. \cr
##'       \code{num.independent.variables} \tab Number of independent variables. \cr
##'       \code{treetype}    \tab Type of forest/tree. Classification, regression or survival. \cr
##'       \code{num.samples}     \tab Number of samples.
##'   }
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @export
predict.ranger <- function(object, data, predict.all = FALSE,
                           seed = NULL, num.threads = NULL,
                           verbose = TRUE, ...) {
  forest <- object$forest
  if (is.null(forest)) {
    stop("Error: No saved forest in ranger object. Please set write.forest to TRUE when calling ranger.")
  }
  predict(forest, data, predict.all, seed, num.threads, verbose)
}
