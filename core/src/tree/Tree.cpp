/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <iterator>
#include "BootstrapSampler.h"

#include "Tree.h"
#include "utility.h"

Tree::Tree(const std::vector<std::vector<size_t>>& child_nodeIDs,
           const std::vector<std::vector<size_t>>& leaf_nodeIDs,
           const std::vector<size_t>& split_varIDs,
           const std::vector<double>& split_values,
           const std::vector<size_t>& oob_sampleIDs,
           const PredictionValues& prediction_values) :
    child_nodeIDs(child_nodeIDs),
    leaf_nodeIDs(leaf_nodeIDs),
    split_varIDs(split_varIDs),
    split_values(split_values),
    oob_sampleIDs(oob_sampleIDs),
    prediction_values(prediction_values) {}

std::vector<size_t> Tree::find_leaf_nodeIDs(Data* prediction_data,
                                            const std::vector<size_t> &sampleIDs) {
  bool use_subsample = !sampleIDs.empty();
  const std::vector<std::vector<size_t>>& child_nodeIDs = get_child_nodeIDs();

  std::vector<size_t> prediction_leaf_nodeIDs;
  prediction_leaf_nodeIDs.resize(prediction_data->get_num_rows());

  size_t num_samples_predict = use_subsample ? sampleIDs.size() : prediction_data->get_num_rows();

  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sampleID = use_subsample ? sampleIDs[i] : i;

    size_t nodeID = 0;
    while (1) {

      // Break if terminal node
      if (child_nodeIDs[0][nodeID] == 0 && child_nodeIDs[1][nodeID] == 0) {
        break;
      }

      // Move to child
      size_t split_varID = get_split_varIDs()[nodeID];
      double value = prediction_data->get(sampleID, split_varID);
      if (value <= get_split_values()[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    }

    prediction_leaf_nodeIDs[sampleID] = nodeID;
  }
  return prediction_leaf_nodeIDs;
}

Tree::~Tree() {}


