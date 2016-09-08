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

##' Case-specific random forests.
##' 
##' In case-specific random forests (CSRF), random forests are built specific to the cases of interest. 
##' Instead of using equal probabilities, the cases are weighted according to their difference to the case of interest.
##' 
##' The algorithm consists of 3 steps: 
##' \enumerate{
##'   \item Grow a random forest on the training data
##'   \item For each observation of interest (test data), the weights of all training observations are computed by counting the number of trees in which both observations are in the same terminal node.
##'   \item For each test observation, grow a weighted random forest on the training data, using the weights obtained in step 2. Predict the outcome of the test observation as usual.
##' }
##'  In total, n+1 random forests are grown, where n is the number observations in the test dataset.
##'  For details, see Xu et al. (2014).
##'
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param training_data Training data of class \code{data.frame}.
##' @param test_data Test data of class \code{data.frame}.
##' @param params1 Parameters for the proximity random forest grown in the first step. 
##' @param params2 Parameters for the prediction random forests grown in the second step. 
##'
##' @return Predictions for the test dataset.
##'
##' @examples
##' ## Split in training and test data
##' train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
##' iris.train <- iris[train.idx, ]
##' iris.test <- iris[-train.idx, ]
##' 
##' ## Run case-specific RF
##' csrf(Species ~ ., training_data = iris.train, test_data = iris.test, 
##'      params1 = list(num.trees = 50, mtry = 4), 
##'      params2 = list(num.trees = 5))
##' 
##' @author Marvin N. Wright
##' @references
##'   Xu, R., Nettleton, D. & Nordman, D.J. (2014). Case-specific random forests. J Comp Graph Stat, in press. DOI: 10.1080/10618600.2014.983641
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



