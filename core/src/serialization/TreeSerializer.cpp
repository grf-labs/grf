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

  Author: Julie Tibshirani (jtibs@cs.stanford.edu)
 #-------------------------------------------------------------------------------*/

#include "Tree.h"
#include "TreeSerializer.h"
#include "utility.h"
#include "PredictionValuesSerializer.h"

void TreeSerializer::serialize(std::ostream& stream, const std::shared_ptr<Tree>& tree) {
  saveVector2D(tree->get_child_nodeIDs(), stream);
  saveVector2D(tree->get_leaf_nodeIDs(), stream);
  saveVector1D(tree->get_split_varIDs(), stream);
  saveVector1D(tree->get_split_values(), stream);
  saveVector1D(tree->get_oob_sampleIDs(), stream);

  PredictionValuesSerializer prediction_values_serializer;
  prediction_values_serializer.serialize(stream, tree->get_prediction_values());
}

std::shared_ptr<Tree> TreeSerializer::deserialize(std::istream& stream) {
  std::vector<std::vector<size_t>> child_nodeIDs;
  readVector2D(child_nodeIDs, stream);

  std::vector<std::vector<size_t>> sampleIDs;
  readVector2D(sampleIDs, stream);

  std::vector<size_t> split_varIDs;
  readVector1D(split_varIDs, stream);

  std::vector<double> split_values;
  readVector1D(split_values, stream);

  std::vector<size_t> oob_sampleIDs;
  readVector1D(oob_sampleIDs, stream);

  PredictionValuesSerializer prediction_values_serializer;
  PredictionValues prediction_values = prediction_values_serializer.deserialize(stream);

  return std::shared_ptr<Tree>(
      new Tree(child_nodeIDs,
               sampleIDs,
               split_varIDs,
               split_values,
               oob_sampleIDs,
               prediction_values));
}
