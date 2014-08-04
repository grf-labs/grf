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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

#include "utility.h"
#include "TreeSurvival.h"
#include "Data.h"

TreeSurvival::TreeSurvival(std::vector<double>* unique_timepoints, size_t status_varID) :
    status_varID(status_varID), unique_timepoints(unique_timepoints) {
  this->num_timepoints = unique_timepoints->size();
}

TreeSurvival::TreeSurvival(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<std::vector<double>> chf, std::vector<double>* unique_timepoints) :
    Tree(child_nodeIDs, split_varIDs, split_values), status_varID(0), unique_timepoints(unique_timepoints), chf(chf) {
  this->num_timepoints = unique_timepoints->size();
}

TreeSurvival::~TreeSurvival() {
}

void TreeSurvival::addPrediction(size_t nodeID, size_t sampleID) {
  predictions[sampleID] = chf[nodeID];
}

void TreeSurvival::appendToFileInternal(std::ofstream& file) {

  // Convert to vector without empty elements and save
  std::vector<size_t> terminal_nodes;
  std::vector<std::vector<double>> chf_vector;
  for (size_t i = 0; i < chf.size(); ++i) {
    if (!chf[i].empty()) {
      terminal_nodes.push_back(i);
      chf_vector.push_back(chf[i]);
    }
  }
  saveVector1D(terminal_nodes, file);
  saveVector2D(chf_vector, file);
}

bool TreeSurvival::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  return findBestSplitLogRank(nodeID, possible_split_varIDs);
}

void TreeSurvival::createEmptyNodeInternal() {
  chf.push_back(std::vector<double>());
}

double TreeSurvival::computePredictionAccuracyInternal() {

  // Compute summed chf for samples
  std::vector<double> sum_chf;
  for (size_t i = 0; i < predictions.size(); ++i) {
    sum_chf.push_back(std::accumulate(predictions[i].begin(), predictions[i].end(), 0));
  }

  // Return concordance index
  return computeConcordanceIndex(data, sum_chf, dependent_varID, status_varID, oob_sampleIDs);
}

bool TreeSurvival::findBestSplitLogRank(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  double best_logrank = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Number of deaths and samples at risk for each timepoint
  size_t* num_deaths = new size_t[num_timepoints];
  size_t* num_samples_at_risk = new size_t[num_timepoints];
  size_t* num_deaths_left_child = new size_t[num_timepoints];
  size_t* num_samples_at_risk_left_child = new size_t[num_timepoints];

  size_t num_unique_death_times = 0;
  computeDeathCounts(num_deaths, num_samples_at_risk, num_unique_death_times, nodeID);

  // Stop early if no split posssible
  if (sampleIDs[nodeID].size() >= 2 * min_node_size) {

    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {

      // Create possible split values
      std::vector<double> possible_split_values;
      data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

      // Try next variable if all equal for this
      if (possible_split_values.size() < 2) {
        continue;
      }

      // For all possible split values
      for (auto& split_value : possible_split_values) {

        // Compute logrank and use split if better than before
        double logrank = computeLogRankTest(nodeID, varID, split_value, num_deaths, num_samples_at_risk,
            num_deaths_left_child, num_samples_at_risk_left_child, num_unique_death_times);
        if (logrank > best_logrank) {
          best_value = split_value;
          best_varID = varID;
          best_logrank = logrank;
        }
      }
    }
  }

  bool result = false;

  // Stop and save CHF if no good split found (this is terminal node).
  if (best_logrank < 0) {
    std::vector<double> chf_temp;
    double chf_value = 0;
    for (size_t i = 0; i < num_timepoints; ++i) {
      if (num_samples_at_risk[i] != 0) {
        chf_value += (double) num_deaths[i] / (double) num_samples_at_risk[i];
      }
      chf_temp.push_back(chf_value);
    }
    chf[nodeID] = chf_temp;
    result = true;
  } else {
    // If not terminal node save best values
    split_varIDs[nodeID] = best_varID;
    split_values[nodeID] = best_value;
  }

  // Clean up
  delete[] num_deaths;
  delete[] num_samples_at_risk;
  delete[] num_deaths_left_child;
  delete[] num_samples_at_risk_left_child;

  return result;
}

void TreeSurvival::computeDeathCounts(size_t* num_deaths, size_t* num_samples_at_risk, size_t& num_unique_death_times,
    size_t nodeID) {

  // Initialize
  for (size_t i = 0; i < num_timepoints; ++i) {
    num_deaths[i] = 0;
    num_samples_at_risk[i] = 0;
  }

  for (auto& sampleID : sampleIDs[nodeID]) {
    double survival_time = data->get(sampleID, dependent_varID);

    size_t t = 0;
    while (t < unique_timepoints->size() && (*unique_timepoints)[t] < survival_time) {
      ++num_samples_at_risk[t];
      ++t;
    }

    // Now t is the survival time, add to at risk and to death if death
    if (t < unique_timepoints->size()) {
      ++num_samples_at_risk[t];
      if (data->get(sampleID, status_varID) == 1) {
        ++num_deaths[t];
      }
    }
  }

  // Count unique death times
  for (size_t j = 0; j < num_timepoints; ++j) {
    if (num_deaths[j] > 0) {
      ++num_unique_death_times;
    }
  }

}

double TreeSurvival::computeLogRankTest(size_t nodeID, size_t varID, double split_value, size_t* num_deaths,
    size_t* num_samples_at_risk, size_t* num_deaths_left_child, size_t* num_samples_at_risk_left_child,
    size_t num_unique_death_times) {

  double nominator = 0;
  double denominator_squared = 0;

  // Initialize
  for (size_t i = 0; i < num_timepoints; ++i) {
    num_deaths_left_child[i] = 0;
    num_samples_at_risk_left_child[i] = 0;
  }

  size_t num_samples_left_child = 0;
  computeChildDeathCounts(nodeID, varID, split_value, num_deaths_left_child, num_samples_at_risk_left_child,
      &num_samples_left_child);

  // Do not consider this split point if fewer than min_node_size samples in one node
  size_t num_samples_right_child = sampleIDs[nodeID].size() - num_samples_left_child;
  if (num_samples_left_child < min_node_size || num_samples_right_child < min_node_size) {
    return -1;
  }

  for (size_t t = 0; t < num_timepoints; ++t) {

    if (num_samples_at_risk[t] < 2) {
      continue;
    }

    // Nominator and demoninator for log-rank test, notation from Ishwaran et al.
    double di = (double) num_deaths[t];
    double di1 = (double) num_deaths_left_child[t];
    double Yi = (double) num_samples_at_risk[t];
    double Yi1 = (double) num_samples_at_risk_left_child[t];
    nominator += di1 - Yi1 * (di / Yi);
    denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;
  }

  if (denominator_squared != 0) {
    return (fabs(nominator / sqrt(denominator_squared)));
  } else {
    return -1;
  }

}

void TreeSurvival::computeChildDeathCounts(size_t nodeID, size_t varID, double split_value,
    size_t* num_deaths_left_child, size_t* num_samples_at_risk_left_child, size_t* num_samples_left_child) {

  // Count deaths and samples at risk in left child for this split
  for (auto& sampleID : sampleIDs[nodeID]) {
    if (data->get(sampleID, varID) <= split_value) {

      ++(*num_samples_left_child);
      double survival_time = data->get(sampleID, dependent_varID);

      size_t t = 0;
      while (t < unique_timepoints->size() && (*unique_timepoints)[t] < survival_time) {
        ++num_samples_at_risk_left_child[t];
        ++t;
      }

      // Now t is the survival time, add to at risk and to death if death
      if (t < unique_timepoints->size()) {
        ++num_samples_at_risk_left_child[t];
        if (data->get(sampleID, status_varID) == 1) {
          ++num_deaths_left_child[t];
        }
      }

    }
  }
}
