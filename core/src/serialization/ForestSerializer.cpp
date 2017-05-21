/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "serialization/ForestSerializer.h"
#include "commons/utility.h"


void ForestSerializer::serialize(std::ostream& stream, const Forest& forest) {
  size_t num_trees = forest.get_trees().size();
  stream.write((char*) &num_trees, sizeof(num_trees));

  for (std::shared_ptr<Tree> tree : forest.get_trees()) {
    tree_serializer.serialize(stream, tree);
  }

  observations_serializer.serialize(stream, forest.get_observations());
}

Forest ForestSerializer::deserialize(std::istream& stream) {
  size_t num_trees;
  stream.read((char*) &num_trees, sizeof(num_trees));
  std::vector<std::shared_ptr<Tree>> trees(num_trees);

  for (size_t i = 0; i < num_trees; ++i) {
    trees[i] = tree_serializer.deserialize(stream);
  }

  Observations observations = observations_serializer.deserialize(stream);
  return Forest(trees, observations);
}

