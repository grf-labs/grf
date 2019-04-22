/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#include <math.h>

#include "SplitFrequencyComputer.h"

std::vector<std::vector<size_t>> SplitFrequencyComputer::compute(const Forest& forest,
                                                                 size_t max_depth) {
  size_t num_variables = forest.get_num_variables();
  std::vector<std::vector<size_t>> result(max_depth, std::vector<size_t>(num_variables));

  for (std::shared_ptr<Tree> tree : forest.get_trees()) {
    const std::vector<std::vector<size_t>>& child_nodes = tree->get_child_nodes();

    size_t depth = 0;
    std::vector<size_t> level = {tree->get_root_node()};

    while (level.size() > 0 && depth < max_depth) {
      std::vector<size_t> next_level;

      for (size_t node : level) {
        if (tree->is_leaf(node)) {
          continue;
        }

        size_t variable = tree->get_split_vars().at(node);
        result[depth][variable]++;

        next_level.push_back(child_nodes[0][node]);
        next_level.push_back(child_nodes[1][node]);
      }

      level = next_level;
      depth++;
    }
  }
  return result;
}

std::map<std::vector<size_t>, size_t> SplitFrequencyComputer::compute_interaction(const Forest& forest) {
    size_t num_variables = forest.get_num_variables();
    std::map<std::vector<size_t>, size_t> result;
    for (std::shared_ptr<Tree> tree : forest.get_trees()) {
        const std::vector<std::vector<size_t>>& child_nodes = tree->get_child_nodes();
        
        std::vector<size_t> level = {tree->get_root_node()};
        
        while (level.size() > 0) {
            std::vector<size_t> next_level;
            
            for (size_t node : level) {
                if (tree->is_leaf(node)) {
                    continue;
                }
                
                std::vector<size_t> nodes_left;
                std::vector<size_t> nodes_right;
    
                size_t variable_root = tree->get_split_vars().at(node);
                size_t variable_left = 0;
                size_t variable_right = 0;
                
                if (tree->is_leaf(child_nodes[0][node])==false) {
                    variable_left = tree->get_split_vars().at(child_nodes[0][node]);
                }
                
                if (tree->is_leaf(child_nodes[1][node])==false) {
                    variable_right = tree->get_split_vars().at(child_nodes[1][node]);
                }
                
                nodes_left.push_back(variable_root);
                nodes_left.push_back(variable_left);
                sort(nodes_left.begin(), nodes_left.end()); //sort the vector
                
                nodes_right.push_back(variable_root);
                nodes_right.push_back(variable_right);
                sort(nodes_right.begin(), nodes_right.end()); //sort the vector
                
                result[nodes_left]++;
                result[nodes_right]++;
                
                next_level.push_back(child_nodes[0][node]);
                next_level.push_back(child_nodes[1][node]);
            }
            
            level = next_level;
        }
    }
    return result;
}
