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
    status_varID(status_varID), unique_timepoints(unique_timepoints), num_deaths(0), num_samples_at_risk(0), num_deaths_left_child(
        0), num_samples_at_risk_left_child(0), num_deaths_0(0), num_samples_at_risk_0(0), num_deaths_1(0), num_samples_at_risk_1(
        0) {
  this->num_timepoints = unique_timepoints->size();
}

TreeSurvival::TreeSurvival(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<std::vector<double>> chf, std::vector<double>* unique_timepoints) :
    Tree(child_nodeIDs, split_varIDs, split_values), status_varID(0), unique_timepoints(unique_timepoints), chf(chf), num_deaths(
        0), num_samples_at_risk(0), num_deaths_left_child(0), num_samples_at_risk_left_child(0), num_deaths_0(0), num_samples_at_risk_0(
        0), num_deaths_1(0), num_samples_at_risk_1(0) {
  this->num_timepoints = unique_timepoints->size();
}

TreeSurvival::~TreeSurvival() {
}

void TreeSurvival::initInternal() {

  // Number of deaths and samples at risk for each timepoint
  num_deaths = new size_t[num_timepoints];
  num_samples_at_risk = new size_t[num_timepoints];
  num_deaths_left_child = new size_t[num_timepoints];
  num_samples_at_risk_left_child = new size_t[num_timepoints];

  // For GWA mode, reuse left for 0
  num_deaths_0 = num_deaths_left_child;
  num_samples_at_risk_0 = num_samples_at_risk_left_child;
  num_deaths_1 = new size_t[num_timepoints];
  num_samples_at_risk_1 = new size_t[num_timepoints];
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
  bool result = false;

  if (splitrule == LOGRANK) {
    result = findBestSplitLogRank(nodeID, possible_split_varIDs);
  } else if (splitrule == AUC || splitrule == AUC_IGNORE_TIES) {
    result = findBestSplitAUC(nodeID, possible_split_varIDs);
  }

  return result;
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

  size_t num_unique_death_times = 0;
  computeDeathCounts(num_unique_death_times, nodeID);

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

      // If the column consists only of gwa data
      if ((possible_split_values.size() == 2 && possible_split_values[0] == 0 && possible_split_values[1] == 1)
          || (possible_split_values.size() == 3 && possible_split_values[0] == 0 && possible_split_values[1] == 1
              && possible_split_values[2] == 2)) {
        findBestSplitValueLogRankGWA(nodeID, varID, best_value, best_varID, best_logrank);
      } else {

        // For all possible split values
        for (auto& split_value : possible_split_values) {

          // Compute logrank and use split if better than before
          double logrank = computeLogRankTest(nodeID, varID, split_value, num_unique_death_times);
          if (logrank > best_logrank) {
            best_value = split_value;
            best_varID = varID;
            best_logrank = logrank;
          }
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

  return result;
}

// TODO: Profile
// TODO: Redundancies with logrank splitting -> fix when merging in master
bool TreeSurvival::findBestSplitAUC(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  double best_auc = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Number of deaths and samples at risk for each timepoint
  size_t* num_deaths = new size_t[num_timepoints];
  size_t* num_samples_at_risk = new size_t[num_timepoints];

  size_t num_unique_death_times = 0;
  computeDeathCounts(num_unique_death_times, nodeID);

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
      double auc = computeAucSplit(nodeID, varID, split_value);

      if (auc > best_auc) {
        best_value = split_value;
        best_varID = varID;
        best_auc = auc;
      }
    }
  }

  bool result = false;

  // Stop and save CHF if no good split found (this is terminal node).
  if (best_auc < 0) {
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

  return result;
}

void TreeSurvival::computeDeathCounts(size_t& num_unique_death_times, size_t nodeID) {

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
      if (data->get(sampleID, status_varID) == 1) {
        ++num_samples_at_risk[t];
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

double TreeSurvival::computeLogRankTest(size_t nodeID, size_t varID, double split_value,
    size_t num_unique_death_times) {

  double nominator = 0;
  double denominator_squared = 0;

  // Initialize
  for (size_t i = 0; i < num_timepoints; ++i) {
    num_deaths_left_child[i] = 0;
    num_samples_at_risk_left_child[i] = 0;
  }

  size_t num_samples_left_child = 0;
  computeChildDeathCounts(nodeID, varID, split_value, &num_samples_left_child);

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
    size_t* num_samples_left_child) {

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

      // Now t is the survival time, add to at risk and to death if death.
      if (t < unique_timepoints->size()) {
        if (data->get(sampleID, status_varID) == 1) {
          ++num_samples_at_risk_left_child[t];
          ++num_deaths_left_child[t];
        }
      }

    }
  }
}

double TreeSurvival::computeAucSplit(size_t nodeID, size_t varID, double split_value) {

  size_t num_node_samples = sampleIDs[nodeID].size();

  // Initialize
  size_t num_total = 0;
  size_t num_left_right = 0;
  size_t num_left_left = 0;
  size_t num_right_right = 0;
  size_t num_samples_left_child = 0;

  for (size_t k = 0; k < num_node_samples; ++k) {
    size_t sample_k = sampleIDs[nodeID][k];
    double time_k = data->get(sample_k, dependent_varID);
    double status_k = data->get(sample_k, status_varID);
    double value_k = data->get(sample_k, varID);

    // Count samples in left node
    if (value_k <= split_value) {
      ++num_samples_left_child;
    }

    for (size_t l = k + 1; l < num_node_samples; ++l) {
      size_t sample_l = sampleIDs[nodeID][l];
      double time_l = data->get(sample_l, dependent_varID);
      double status_l = data->get(sample_l, status_varID);
      double value_l = data->get(sample_l, varID);

      // Assign smaller and larger values, do not count if times equal
      double value_smaller;
      double value_larger;
      double status_smaller;

      if (time_k < time_l) {
        value_smaller = value_k;
        value_larger = value_l;
        status_smaller = status_k;
      } else if (time_l < time_k) {
        value_smaller = value_l;
        value_larger = value_k;
        status_smaller = status_l;
      } else {
        if (splitrule != AUC_IGNORE_TIES) {
          if (time_k == time_l && status_k == 1 && status_l == 1 && value_k != value_l) {
            // Tie in survival time, event for both and no tie for covariate -> count with 0.5
            ++num_total;
            ++num_left_left;
          }
        }
        continue;
      }

      // Do not count if smaller time censored
      if (status_smaller == 0) {
        continue;
      }

      // Count total
      ++num_total;

      // Count
      if (value_smaller <= split_value && value_larger > split_value) {
        ++num_left_right;
      } else if (value_smaller <= split_value && value_larger <= split_value) {
        ++num_left_left;
      } else if (value_smaller > split_value && value_larger > split_value) {
        ++num_right_right;
      }

    }
  }

  // Do not consider this split point if fewer than min_node_size samples in one node
  size_t num_samples_right_child = num_node_samples - num_samples_left_child;
  if (num_samples_left_child < min_node_size || num_samples_right_child < min_node_size) {
    return -1;
  } else {
    return fabs((num_left_right + 0.5 * num_left_left + 0.5 * num_right_right) / num_total - 0.5);
  }
}

void TreeSurvival::findBestSplitValueLogRankGWA(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
    double& best_logrank) {

  double nominator_split_0 = 0;
  double nominator_split_1 = 0;
  double denominator_squared_split_0 = 0;
  double denominator_squared_split_1 = 0;
  size_t all_samples_0 = 0;
  size_t all_samples_1 = 0;

  size_t* num_samples_at_risk_value = num_samples_at_risk_0;
  size_t* num_deaths_value = num_deaths_0;

  // Initialize
  for (size_t i = 0; i < num_timepoints; ++i) {
    num_samples_at_risk_0[i] = 0;
    num_samples_at_risk_1[i] = 0;

    num_deaths_0[i] = 0;
    num_deaths_1[i] = 0;
  }

  // Count deaths and samples at risk in left child for this split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    double survival_time = data->get(sampleID, dependent_varID);

    if (value == 0) {
      ++all_samples_0;
      num_samples_at_risk_value = num_samples_at_risk_0;
      num_deaths_value = num_deaths_0;
    } else if (value == 1) {
      ++all_samples_1;
      num_samples_at_risk_value = num_samples_at_risk_1;
      num_deaths_value = num_deaths_1;
    }

    size_t t = 0;
    while (t < unique_timepoints->size() && (*unique_timepoints)[t] < survival_time) {
      ++num_samples_at_risk_value[t];
      ++t;
    }

    // Now t is the survival time, add to at risk and to death if death
    if (t < unique_timepoints->size()) {
      if (data->get(sampleID, status_varID) == 1) {
        ++num_samples_at_risk_value[t];
        ++num_deaths_value[t];
      }
    }
  }

  // Do not consider this split point if fewer than min_node_size samples in one node
  size_t all_samples = sampleIDs[nodeID].size();
  bool split_okay_0 = true;
  if (all_samples_0 < min_node_size || (all_samples - all_samples_0) < min_node_size) {
    split_okay_0 = false;
  }
  bool split_okay_1 = true;
  if ((all_samples_0 + all_samples_1) < min_node_size
      || (all_samples - all_samples_0 - all_samples_1) < min_node_size) {
    split_okay_1 = false;
  }
  if (!split_okay_0 && !split_okay_1) {
    return;
  }

  for (size_t t = 0; t < num_timepoints; ++t) {

    if (num_samples_at_risk[t] < 2) {
      continue;
    }

    // Nominator and demoninator for log-rank test, notation from Ishwaran et al.
    double di = (double) num_deaths[t];
    double Yi = (double) num_samples_at_risk[t];
    double temp = ((Yi - di) / (Yi - 1)) * di;

    if (split_okay_0) {
      double di1_split_0 = (double) num_deaths_0[t];
      double Yi1_split_0 = (double) num_samples_at_risk_0[t];
      nominator_split_0 += di1_split_0 - Yi1_split_0 * (di / Yi);
      denominator_squared_split_0 += (Yi1_split_0 / Yi) * (1.0 - Yi1_split_0 / Yi) * temp;
    }

    if (split_okay_1) {
      double di1_split_1 = (double) (num_deaths_0[t] + num_deaths_1[t]);
      double Yi1_split_1 = (double) (num_samples_at_risk_0[t] + num_samples_at_risk_1[t]);
      nominator_split_1 += di1_split_1 - Yi1_split_1 * (di / Yi);
      denominator_squared_split_1 += (Yi1_split_1 / Yi) * (1.0 - Yi1_split_1 / Yi) * temp;
    }
  }
}

