/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "tree/Tree.h"
#include "serialization/TreeSerializer.h"
#include "commons/utility.h"
#include "serialization/PredictionValuesSerializer.h"

void TreeSerializer::serialize(std::ostream& stream, const std::shared_ptr<Tree>& tree) {
  size_t root_node = tree->get_root_node();
  stream.write((char*) &root_node, sizeof(root_node));

  write_matrix(tree->get_child_nodes(), stream);
  write_matrix(tree->get_leaf_samples(), stream);
  write_vector(tree->get_split_vars(), stream);
  write_vector(tree->get_split_values(), stream);
  write_vector(tree->get_oob_samples(), stream);

  PredictionValuesSerializer prediction_values_serializer;
  prediction_values_serializer.serialize(stream, tree->get_prediction_values());
}

std::shared_ptr<Tree> TreeSerializer::deserialize(std::istream& stream) {
  size_t root_node;
  stream.read((char*) &root_node, sizeof(root_node));

  std::vector<std::vector<size_t>> child_nodes;
  read_matrix(child_nodes, stream);

  std::vector<std::vector<size_t>> samples;
  read_matrix(samples, stream);

  std::vector<size_t> split_vars;
  read_vector(split_vars, stream);

  std::vector<double> split_values;
  read_vector(split_values, stream);

  std::vector<size_t> oob_samples;
  read_vector(oob_samples, stream);

  PredictionValuesSerializer prediction_values_serializer;
  PredictionValues prediction_values = prediction_values_serializer.deserialize(stream);

  return std::shared_ptr<Tree>(
      new Tree(root_node,
               child_nodes,
               samples,
               split_vars,
               split_values,
               oob_samples,
               prediction_values));
}
