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

TreeSurvival::TreeSurvival(std::vector<double>* unique_timepoints, size_t status_varID,
    std::vector<size_t>* response_timepointIDs) :
    status_varID(status_varID), unique_timepoints(unique_timepoints), response_timepointIDs(response_timepointIDs), num_deaths(
        0), num_samples_at_risk(0) {
  this->num_timepoints = unique_timepoints->size();
}

TreeSurvival::TreeSurvival(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<std::vector<double>> chf, std::vector<double>* unique_timepoints,
    std::vector<size_t>* response_timepointIDs, std::vector<bool>* is_ordered_variable) :
    Tree(child_nodeIDs, split_varIDs, split_values, is_ordered_variable), status_varID(0), unique_timepoints(
        unique_timepoints), response_timepointIDs(response_timepointIDs), chf(chf), num_deaths(0), num_samples_at_risk(
        0) {
  this->num_timepoints = unique_timepoints->size();
}

TreeSurvival::~TreeSurvival() {
}

void TreeSurvival::initInternal() {
  // Number of deaths and samples at risk for each timepoint
  num_deaths = new size_t[num_timepoints];
  num_samples_at_risk = new size_t[num_timepoints];
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
  for (size_t i = 0; i < prediction_terminal_nodeIDs.size(); ++i) {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[i];
    sum_chf.push_back(std::accumulate(chf[terminal_nodeID].begin(), chf[terminal_nodeID].end(), 0));
  }

  // Return concordance index
  return computeConcordanceIndex(data, sum_chf, dependent_varID, status_varID, oob_sampleIDs);
}

bool TreeSurvival::findBestSplitLogRank(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = sampleIDs[nodeID].size();
  double best_logrank = -1;
  size_t best_varID = 0;
  double best_value = 0;

  computeDeathCounts(nodeID);

  // Stop early if no split posssible
  if (num_samples_node >= 2 * min_node_size) {

    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {

      // Find best split value, if ordered consider all values as split values, else all 2-partitions
      if ((*is_ordered_variable)[varID]) {

        // TODO: Here

        double best_value_temp = best_value;
        size_t best_varID_temp = best_varID;
        double best_logrank_temp = best_logrank;

//        std::cout << "Old: " << std::endl;
        findBestSplitValueLogRank(nodeID, varID, best_value, best_varID, best_logrank);

//        std::cout << "New: " << std::endl;
        findBestSplitValueLogRankNew(nodeID, varID, num_samples_node, best_value_temp, best_varID_temp,
            best_logrank_temp);

//        std::cout << std::endl;
//        std::getchar();

        if (best_logrank != best_logrank_temp && best_value != best_value_temp) {
          std::cout << "best_value: " << best_value << ", " << best_value_temp << std::endl;
          std::cout << "best_varID: " << best_varID << ", " << best_varID_temp << std::endl;
          std::cout << "best_logrank: " << best_logrank << ", " << best_logrank_temp << std::endl;
          std::cout << std::endl;
        }

//        findBestSplitValueLogRank(nodeID, varID, best_value, best_varID, best_logrank);
//        findBestSplitValueLogRankNew(nodeID, varID, num_samples_node, best_value, best_varID, best_logrank);

      } else {
        findBestSplitValueLogRankUnordered(nodeID, varID, best_value, best_varID, best_logrank);
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

void TreeSurvival::computeDeathCounts(size_t nodeID) {

  // Initialize
  for (size_t i = 0; i < num_timepoints; ++i) {
    num_deaths[i] = 0;
    num_samples_at_risk[i] = 0;
  }

  for (auto& sampleID : sampleIDs[nodeID]) {
    double survival_time = data->get(sampleID, dependent_varID);

    size_t t = 0;
    while (t < num_timepoints && (*unique_timepoints)[t] < survival_time) {
      ++num_samples_at_risk[t];
      ++t;
    }

    // Now t is the survival time, add to at risk and to death if death
    if (t < num_timepoints) {
      if (data->get(sampleID, status_varID) == 1) {
        ++num_samples_at_risk[t];
        ++num_deaths[t];
      }
    }
  }
}

void TreeSurvival::computeChildDeathCounts(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
    size_t* num_samples_right_child, size_t* delta_samples_at_risk_right_child, size_t* num_deaths_right_child) {
  size_t num_splits = possible_split_values.size();

  // Count deaths in right child per timepoint and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    size_t survival_timeID = (*response_timepointIDs)[sampleID];

    // Count deaths until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {

      if (value > possible_split_values[i]) {
        ++num_samples_right_child[i];
        ++delta_samples_at_risk_right_child[i * num_timepoints + survival_timeID];
        if (data->get(sampleID, status_varID) == 1) {
          ++num_deaths_right_child[i * num_timepoints + survival_timeID];
        }
      } else {
        break;
      }
    }
  }
}

void TreeSurvival::findBestSplitValueLogRank(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
    double& best_logrank) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // Remove largest value because no split possible
  possible_split_values.pop_back();

  size_t num_splits = possible_split_values.size();

  // Initialize
  size_t* num_deaths_right_child = new size_t[num_splits * num_timepoints]();
  size_t* delta_samples_at_risk_right_child = new size_t[num_splits * num_timepoints]();
  size_t* num_samples_right_child = new size_t[num_splits]();

  computeChildDeathCounts(nodeID, varID, possible_split_values, num_samples_right_child,
      delta_samples_at_risk_right_child, num_deaths_right_child);

//  // TODO: Remove
//  std::cout << "delta_at_risk old: " << std::endl;
//  for (size_t i = 0; i < num_splits; ++i) {
//    std::cout << possible_split_values[i] << ": ";
//    for (size_t t = 0; t < num_timepoints; ++t) {
//      std::cout << delta_samples_at_risk_right_child[i * num_timepoints + t] << " ";
//    }
//    std::cout << std::endl;
//  }
//  std::cout << std::endl << std::endl;

//// TODO: Remove
//  std::cout << "death old: " << std::endl;
//  for (size_t i = 0; i < num_splits; ++i) {
//    std::cout << possible_split_values[i] << ": ";
//    for (size_t t = 0; t < num_timepoints; ++t) {
//      std::cout << num_deaths_right_child[i * num_timepoints + t] << " ";
//    }
//    std::cout << std::endl;
//  }
//  std::cout << std::endl << std::endl;

// Compute logrank test for all splits and use best
  for (size_t i = 0; i < num_splits; ++i) {
    double numerator = 0;
    double denominator_squared = 0;

    // Stop if minimal node size reached
    size_t num_samples_left_child = sampleIDs[nodeID].size() - num_samples_right_child[i];
    if (num_samples_right_child[i] < min_node_size || num_samples_left_child < min_node_size) {
      continue;
    }

//    // TODO: Remove
//    std::cout << "split: " << possible_split_values[i] << std::endl;
//    std::cout << "old at risk ";

// Compute logrank test statistic for this split
    size_t num_samples_at_risk_right_child = num_samples_right_child[i];
    for (size_t t = 0; t < num_timepoints; ++t) {
      if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
        break;
      }

      if (num_deaths[t] > 0) {
        // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
        double di = (double) num_deaths[t];
        double di1 = (double) num_deaths_right_child[i * num_timepoints + t];
        double Yi = (double) num_samples_at_risk[t];
        double Yi1 = (double) num_samples_at_risk_right_child;
        numerator += di1 - Yi1 * (di / Yi);
        denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;

        //TODO: Remove
//        std::cout << "time " << t << ", di " << di << ", di1 " << di1 << ", Yi " << Yi << ", Yi1 " << Yi1
//            << ", numerator " << numerator << ", denominator_squared " << denominator_squared << std::endl;
      }

//      // TODO: Remove
//      std::cout << num_samples_at_risk_right_child << " ";

// Reduce number of samples at risk for next timepoint
      num_samples_at_risk_right_child -= delta_samples_at_risk_right_child[i * num_timepoints + t];

    }
    double logrank = -1;
    if (denominator_squared != 0) {
      logrank = fabs(numerator / sqrt(denominator_squared));
    }

//    //TODO: Remove
//    std::cout << "split: " << possible_split_values[i] << ", logrank " << logrank << ", numerator " << numerator
//        << ", denominator_squared " << denominator_squared << std::endl;

    if (logrank > best_logrank) {
      best_value = possible_split_values[i];
      best_varID = varID;
      best_logrank = logrank;
    }
  }

  delete[] num_deaths_right_child;
  delete[] delta_samples_at_risk_right_child;
  delete[] num_samples_right_child;
}

// TODO: Try without extra loop ..
// TODO: Use left counts instead of right?
void TreeSurvival::findBestSplitValueLogRankNew(size_t nodeID, size_t varID, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_logrank) {

  // TODO: Cleanup
  // TODO: Stop if only one unique value?
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::vector<size_t> count(num_unique);
  std::vector<size_t> death(num_unique * num_timepoints);
  std::vector<size_t> delta_at_risk(num_unique * num_timepoints);

  for (auto& sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);
    size_t survival_timeID = (*response_timepointIDs)[sampleID];

    ++count[index];
    ++delta_at_risk[index * num_timepoints + survival_timeID];
    if (data->get(sampleID, status_varID) == 1) {
      ++death[index * num_timepoints + survival_timeID];
    }
  }

  // TODO: Remove
//  std::cout << "delta_at_risk new: " << std::endl;
//  for (size_t i = 0; i < num_unique; ++i) {
//    for (size_t t = 0; t < num_timepoints; ++t) {
//      std::cout << delta_at_risk[i * num_timepoints + t] << " ";
//    }
//    std::cout << std::endl;
//  }
//  std::cout << std::endl << std::endl;

  // TODO: Dont use extra loop?
  std::vector<size_t> at_risk((num_unique - 1) * num_timepoints);
  std::vector<size_t> death2(num_unique * num_timepoints);
  for (size_t t = 0; t < num_timepoints; ++t) {
    size_t sum_at_risk = 0;
    size_t sum_death = 0;
    for (size_t i = 0; i < num_unique - 1; ++i) {
      size_t ii = num_unique - 2 - i;
      sum_at_risk += delta_at_risk[(ii + 1) * num_timepoints + t];

//      std::cout << sum_at_risk << " ";

      at_risk[ii * num_timepoints + t] = sum_at_risk;
      sum_death += death[(ii + 1) * num_timepoints + t];
      death2[ii * num_timepoints + t] = sum_death;
    }
//    std::cout << std::endl;
  }

  // TODO: -1?
  std::vector<size_t> death3((num_unique-1) * num_timepoints);
  for (size_t t = 0; t < num_timepoints; ++t) {
    size_t sum_death = 0;
    for (size_t i = 0; i < num_unique-1; ++i) {
      sum_death += death[i * num_timepoints + t];
      death3[i * num_timepoints + t] = sum_death;
      if (death3[i * num_timepoints + t] != num_deaths[t] - death2[i * num_timepoints + t]) {
        std::cout << "DIFF death" << std::endl;
      }
    }
  }

  std::cout << "AAAAAA" << std::endl << std::endl;

  // TODO: Missing: Number of samples at risk in total
  std::vector<size_t> at_risk3((num_unique-1) * num_timepoints);
  for (size_t t = 0; t < num_timepoints; ++t) {
    size_t sum_at_risk = 0;
    for (size_t i = 0; i < num_unique-1; ++i) {

      sum_at_risk += delta_at_risk[i * num_timepoints + t];

      std::cout << num_samples_at_risk[t]-sum_at_risk << " ";

      at_risk3[i * num_timepoints + t] = sum_at_risk;
      //if (at_risk3[i * num_timepoints + t] != num_samples_at_risk[t] - at_risk[i * num_timepoints + t]) {
     //   std::cout << "DIFF at_risk: " << num_samples_at_risk[t] - at_risk[i * num_timepoints + t] << " " << at_risk3[i * num_timepoints + t]  << std::endl;
      //}

    }
    std::cout << std::endl;
  }

  std::cout << "BBBBB" << std::endl << std::endl;

//  // TODO: Remove
//  std::cout << "at_risk new: " << std::endl;
//  for (size_t i = 0; i < num_unique-1; ++i) {
//    std::cout << data->getUniqueDataValue(varID, i) << ": ";
//    for (size_t t = 0; t < num_timepoints; ++t) {
//      std::cout << at_risk[i * num_timepoints + t] << " ";
//    }
//    std::cout << std::endl;
//  }
//  std::cout << std::endl << std::endl;

//  std::cout << "death new: " << std::endl;
//  for (size_t i = 0; i < num_unique - 1; ++i) {
//    std::cout << data->getUniqueDataValue(varID, i) << ": ";
//    for (size_t t = 0; t < num_timepoints; ++t) {
//      std::cout << death2[i * num_timepoints + t] << " ";
//    }
//    std::cout << std::endl;
//  }
//  std::cout << std::endl << std::endl;

  size_t n_left = 0;

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique; ++i) {

    // Stop if nothing here
//    if (count[i] == 0) {
//      continue;
//    }

    n_left += count[i];

    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }

    // Stop if minimal node size reached
    if (n_left < min_node_size || n_right < min_node_size) {
      continue;
    }

    // Init logrank test
    double numerator = 0;
    double denominator_squared = 0;
//    size_t at_risk = n_left;

//    std::cout << "new at risk ";

    // Compute logrank test
    size_t num_samples_at_risk_right_child = n_right;
    for (size_t t = 0; t < num_timepoints; ++t) {

      if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
        break;
      }

//      if (num_samples_at_risk[t] < 2) {
//        continue;
//      }

      // TODO: Dont use minus -> left child
      // TODO: at_risk now counts right child
      if (num_samples_at_risk[t] > 1 && num_deaths[t] > 0) {
        // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
        double di = (double) num_deaths[t];
        //double di1 = di - (double) death2[i * num_timepoints + t];
        double di1 = (double) death3[i * num_timepoints + t];
        double Yi = (double) num_samples_at_risk[t];
        double Yi1 = Yi - (double) num_samples_at_risk_right_child;
        numerator += di1 - Yi1 * (di / Yi);
        denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;

        //TODO: Remove
//        std::cout << "time " << t << ", di " << di << ", di1 " << di1 << ", Yi " << Yi << ", Yi1 " << Yi1
//            << ", numerator " << numerator << ", denominator_squared " << denominator_squared << std::endl;

      }

      // Reduce number of samples at risk for next timepoint
      //at_risk -= delta_at_risk[i * num_timepoints + t];
      num_samples_at_risk_right_child -= at_risk[i * num_timepoints + t];

      // TODO: Remove
      std::cout << num_samples_at_risk_right_child << " ";

    }

    // TODO: Remove
    std::cout << std::endl;

    double logrank = -1;
    if (denominator_squared != 0) {
      logrank = fabs(numerator / sqrt(denominator_squared));
    }

//    //TODO: Remove
//    std::cout << "split: " << data->getUniqueDataValue(varID, i) << ", logrank " << logrank << ", numerator "
//        << numerator << ", denominator_squared " << denominator_squared << std::endl;

    // If better than before, use this
    if (logrank > best_logrank) {
      best_value = data->getUniqueDataValue(varID, i);
      best_varID = varID;
      best_logrank = logrank;
    }
  }
}

void TreeSurvival::findBestSplitValueLogRankUnordered(size_t nodeID, size_t varID, double& best_value,
    size_t& best_varID, double& best_logrank) {

  // Create possible split values
  std::vector<double> factor_levels;
  data->getAllValues(factor_levels, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (factor_levels.size() < 2) {
    return;
  }

  // Number of possible splits is 2^num_levels
  size_t num_splits = (1 << factor_levels.size());

  // Compute logrank test statistic for each possible split
  // Split where all left (0) or all right (1) are excluded
  // The second half of numbers is just left/right switched the first half -> Exclude second half
  for (size_t local_splitID = 1; local_splitID < num_splits / 2; ++local_splitID) {

    // Compute overall splitID by shifting local factorIDs to global positions
    size_t splitID = 0;
    for (size_t j = 0; j < factor_levels.size(); ++j) {
      if ((local_splitID & (1 << j))) {
        double level = factor_levels[j];
        size_t factorID = floor(level) - 1;
        splitID = splitID | (1 << factorID);
      }
    }

    // Initialize
    size_t* num_deaths_right_child = new size_t[num_timepoints]();
    size_t* delta_samples_at_risk_right_child = new size_t[num_timepoints]();
    size_t num_samples_right_child = 0;
    double numerator = 0;
    double denominator_squared = 0;

    // Count deaths in right child per timepoint
    for (auto& sampleID : sampleIDs[nodeID]) {
      size_t survival_timeID = (*response_timepointIDs)[sampleID];
      double value = data->get(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1 << factorID))) {
        ++num_samples_right_child;
        ++delta_samples_at_risk_right_child[survival_timeID];
        if (data->get(sampleID, status_varID) == 1) {
          ++num_deaths_right_child[survival_timeID];
        }
      }

    }

    // Stop if minimal node size reached
    size_t num_samples_left_child = sampleIDs[nodeID].size() - num_samples_right_child;
    if (num_samples_right_child < min_node_size || num_samples_left_child < min_node_size) {
      delete[] num_deaths_right_child;
      delete[] delta_samples_at_risk_right_child;
      continue;
    }

    // Compute logrank test statistic for this split
    size_t num_samples_at_risk_right_child = num_samples_right_child;
    for (size_t t = 0; t < num_timepoints; ++t) {
      if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
        break;
      }

      if (num_deaths[t] > 0) {
        // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
        double di = (double) num_deaths[t];
        double di1 = (double) num_deaths_right_child[t];
        double Yi = (double) num_samples_at_risk[t];
        double Yi1 = (double) num_samples_at_risk_right_child;
        numerator += di1 - Yi1 * (di / Yi);
        denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;
      }

      // Reduce number of samples at risk for next timepoint
      num_samples_at_risk_right_child -= delta_samples_at_risk_right_child[t];
    }
    double logrank = -1;
    if (denominator_squared != 0) {
      logrank = fabs(numerator / sqrt(denominator_squared));
    }

    if (logrank > best_logrank) {
      best_value = splitID;
      best_varID = varID;
      best_logrank = logrank;
    }

    delete[] num_deaths_right_child;
    delete[] delta_samples_at_risk_right_child;
  }

}
