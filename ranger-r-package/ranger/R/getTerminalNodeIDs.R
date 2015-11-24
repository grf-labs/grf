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

##' Get terminal node IDs of observations.
##'
##' @param rf \code{ranger} object.
##' @param dat New dataset. Terminal node IDs for this dataset are obtained. 
##'
##' @return Matrix with terminal nodeIDs for all observations in dataset and trees.
##'
##' @examples
##' library(ranger)
##' rf <- ranger(Species ~ ., data = iris, write.forest = TRUE)
##' getTerminalNodeIDs(rf, iris)
##' @export
getTerminalNodeIDs <- function(rf, dat) {
  ## Check if forests exists
  if (is.null(rf$forest)) {
    stop("Error: No saved forest in ranger object. Please set write.forest to TRUE when calling ranger.")
  }
  
  ## Get terminal node IDs for each observation
  result <- t(apply(dat, 1, function(obs) {
    ## Drop observation down the trees
    sapply(1:rf$num.trees, function(tree) {
      
      child.nodeIDs <- rf$forest$child.nodeIDs[[tree]]
      split.varIDs <- rf$forest$split.varIDs[[tree]]
      split.values <- rf$forest$split.values[[tree]]
      
      nodeID <- 1
      while (1) {
        
        ## Break if terminal node
        if (length(child.nodeIDs[[nodeID]]) < 1) {
          break
        }
        
        ## Move to child
        split.varID <- split.varIDs[[nodeID]]
        value <- obs[split.varID]
        if (value <= split.values[[nodeID]]) {
          ## Move to left child
          nodeID <- child.nodeIDs[[nodeID]][1] + 1
        } else {
          ## Move to right child
          nodeID <- child.nodeIDs[[nodeID]][2] + 1
        }
      }
      
      return(nodeID)
    })
  }))
  
  return(result)
}