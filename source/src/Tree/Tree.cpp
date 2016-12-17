/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 wright@imbs.uni-luebeck.de
 #-------------------------------------------------------------------------------*/

#include <iterator>
#include <Forest/BootstrapSampler.h>

#include "Tree.h"
#include "utility.h"

Tree::Tree(std::vector<std::vector<size_t>> &child_nodeIDs,
           std::vector<std::vector<size_t>> sampleIDs,
           std::vector<size_t> &split_varIDs,
           std::vector<double> &split_values,
           Data *data,
           BootstrapSampler* bootstrap_sampler) :
    child_nodeIDs(child_nodeIDs),
    sampleIDs(sampleIDs),
    split_varIDs(split_varIDs),
    split_values(split_values),
    data(0),
    bootstrap_sampler(bootstrap_sampler) {}

Tree::~Tree() {}

void Tree::predict(const Data* prediction_data, bool oob_prediction) {

  size_t num_samples_predict;
  if (oob_prediction) {
    num_samples_predict = bootstrap_sampler->getNumSamplesOob();
  } else {
    num_samples_predict = prediction_data->getNumRows();
  }

  prediction_terminal_nodeIDs.resize(num_samples_predict, 0);

// For each sample start in root, drop down the tree and return final value
  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sample_idx;
    if (oob_prediction) {
      sample_idx = bootstrap_sampler->getOobSampleIDs()[i];
    } else {
      sample_idx = i;
    }
    size_t nodeID = 0;
    while (1) {

      // Break if terminal node
      if (child_nodeIDs[0][nodeID] == 0 && child_nodeIDs[1][nodeID] == 0) {
        break;
      }

      // Move to child
      size_t split_varID = split_varIDs[nodeID];
      double value = prediction_data->get(sample_idx, split_varID);
      if (value <= split_values[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    }

    prediction_terminal_nodeIDs[i] = nodeID;
  }
}

std::vector<size_t> Tree::get_neighboring_samples(size_t sampleID) {
  size_t nodeID = prediction_terminal_nodeIDs[sampleID];
  return sampleIDs[nodeID];
}

