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

##' Case specific random forests.
##'
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param training_data Training data of class \code{data.frame}.
##' @param test_data Test data of class \code{data.frame}.
##' @param params1 Parameters for the proximity random forest grown in the first step. 
##' @param params2 Parameters for the prediction random forests grown in the second step. 
##'
##' @return Predictions for the test dataset
##'
##' @examples
##' #TODO
##' @export
csrf <- function(formula, training_data, test_data, params1 = list(), params2 = list()) {
  ## Grow a random forest on the training data to obtain weights
  rf.proximity <- do.call(ranger, c(list(formula = formula, data = training_data, 
                                         write.forest = TRUE), params1))
  
  ## Get terminal nodes
  terminal.nodeIDs.train <- getTerminalNodeIDs(rf.proximity, training_data)
  terminal.nodeIDs.test <- getTerminalNodeIDs(rf.proximity, test_data)
  
  ## Grow weighted RFs for test observations, predict the outcome
  predictions <- sapply(1:nrow(test_data), function(i) {
    ## Compute weights from first RF
    num.same.node <- rowSums(terminal.nodeIDs.test[i, ] == terminal.nodeIDs.train)
    weights <- num.same.node / sum(num.same.node)
    
    ## Grow weighted RF
    rf.prediction <- do.call(ranger, c(list(formula = formula, data = training_data, 
                                            write.forest = TRUE, case.weights = weights), 
                                       params2))
    
    ## Predict outcome
    predict(rf.prediction, test_data[i, ])$predictions
  })
  
  ## Return predictions
  predictions
}



