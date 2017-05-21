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

#ifndef GRF_TREE_H_
#define GRF_TREE_H_

#include <vector>
#include <random>
#include <iostream>

#include "commons/globals.h"
#include "commons/Data.h"
#include "sampling/BootstrapSampler.h"
#include "prediction/PredictionValues.h"
#include "splitting/SplittingRule.h"

class Tree {
public:
  Tree(const std::vector<std::vector<size_t>>& child_nodeIDs,
       const std::vector<std::vector<size_t>>& leaf_nodeIDs,
       const std::vector<size_t>& split_varIDs,
       const std::vector<double>& split_values,
       const std::vector<size_t>& oob_sampleIDs,
       const PredictionValues& prediction_values);

  ~Tree();

  std::vector<size_t> find_leaf_nodeIDs(Data *prediction_data,
                                        const std::vector<size_t> &sampleIDs);

  const std::vector<std::vector<size_t>>& get_child_nodeIDs() {
    return child_nodeIDs;
  }

  const std::vector<std::vector<size_t>>& get_leaf_nodeIDs() {
    return leaf_nodeIDs;
  }

  const std::vector<double>& get_split_values() {
    return split_values;
  }

  const std::vector<size_t>& get_split_varIDs() {
    return split_varIDs;
  }

  const std::vector<size_t>& get_oob_sampleIDs() {
    return oob_sampleIDs;
  }

  const PredictionValues& get_prediction_values() {
    return prediction_values;
  }

  void set_leaf_nodeIDs(const std::vector<std::vector<size_t>>& leaf_nodeIDs) {
    this->leaf_nodeIDs = leaf_nodeIDs;
  }

  void set_oob_sampleIDs(std::vector<size_t>& oob_sampleIDs) {
    this->oob_sampleIDs = oob_sampleIDs;
  }

  void set_prediction_values(const PredictionValues& prediction_values) {
    this->prediction_values = prediction_values;
  }

private:
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> leaf_nodeIDs;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  std::vector<size_t> oob_sampleIDs;

  PredictionValues prediction_values;
};

#endif /* GRF_TREE_H_ */
