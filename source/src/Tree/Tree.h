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

#ifndef TREE_H_
#define TREE_H_

#include <vector>
#include <random>
#include <iostream>

#include "globals.h"
#include "Data.h"
#include "BootstrapSampler.h"
#include "SplittingRule.h"
#include "RelabelingStrategy.h"

class Tree {
public:
  Tree(std::vector<std::vector<size_t>>& child_nodeIDs,
       std::vector<std::vector<size_t>> sampleIDs,
       std::vector<size_t>& split_varIDs,
       std::vector<double>& split_values,
       Data* data,
       BootstrapSampler* bootstrap_sampler);

  ~Tree();

  void predict(const Data* prediction_data, bool oob_prediction);

  const std::vector<std::vector<size_t> >& getChildNodeIDs() const {
    return child_nodeIDs;
  }

  const std::vector<std::vector<size_t>>& get_sampleIDs() const {
    return sampleIDs;
  }

  const std::vector<double>& getSplitValues() const {
    return split_values;
  }

  const std::vector<size_t>& getSplitVarIDs() const {
    return split_varIDs;
  }

  const std::vector<size_t>& getOobSampleIDs() const {
    return bootstrap_sampler->getOobSampleIDs();
  }

  const std::vector<size_t>& getInbagCounts() const {
    return bootstrap_sampler->getInbagCounts();
  }

  std::vector<size_t> get_neighboring_samples(size_t sampleID);

private:
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> sampleIDs;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  Data* data;

  BootstrapSampler* bootstrap_sampler;

  // When growing here the OOB set is used
  // Terminal nodeIDs for prediction samples
  std::vector<size_t> prediction_terminal_nodeIDs;

private:
  DISALLOW_COPY_AND_ASSIGN(Tree);
};

#endif /* TREE_H_ */
