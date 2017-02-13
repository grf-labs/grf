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

#include "TreePredictor.h"

std::vector<size_t> TreePredictor::get_terminal_nodeIDs(std::shared_ptr<Tree> tree,
                                                        Data* prediction_data,
                                                        const std::vector<size_t>& sampleIDs) {
  bool use_subsample = !sampleIDs.empty();
  const std::vector<std::vector<size_t>>& child_nodeIDs = tree->get_child_nodeIDs();
  const std::vector<double>& split_values = tree->get_split_values();
  const std::vector<size_t>& split_varIDs = tree->get_split_varIDs();

  std::vector<size_t> prediction_terminal_nodeIDs;
  prediction_terminal_nodeIDs.resize(prediction_data->get_num_rows());

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
      size_t split_varID = split_varIDs[nodeID];
      double value = prediction_data->get(sampleID, split_varID);
      if (value <= split_values[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    }

    prediction_terminal_nodeIDs[sampleID] = nodeID;
  }
  return prediction_terminal_nodeIDs;
}
