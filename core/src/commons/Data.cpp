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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <iterator>
#include <stdexcept>

#include "Data.h"

namespace grf {

Data::Data(const double* data_ptr, size_t num_rows, size_t num_cols) {
  if (data_ptr == nullptr) {
    throw std::runtime_error("Invalid data storage: nullptr");
  }
  this->data_ptr = data_ptr;
  this->num_rows = num_rows;
  this->num_cols = num_cols;
}

Data::Data(const std::vector<double>& data, size_t num_rows, size_t num_cols) :
  Data(data.data(), num_rows, num_cols) {}

Data::Data(const std::pair<std::vector<double>, std::vector<size_t>>& data) :
  Data(data.first.data(), data.second.at(0), data.second.at(1)) {}

void Data::set_outcome_index(size_t index) {
  set_outcome_index(std::vector<size_t>({index}));
}

void Data::set_outcome_index(const std::vector<size_t>& index) {
  this->outcome_index = index;
  disallowed_split_variables.insert(index.begin(), index.end());
}

void Data::set_treatment_index(size_t index) {
  set_treatment_index(std::vector<size_t>({index}));
}

void Data::set_treatment_index(const std::vector<size_t>& index) {
  this->treatment_index = index;
  disallowed_split_variables.insert(index.begin(), index.end());
}

void Data::set_instrument_index(size_t index) {
  this->instrument_index = index;
  disallowed_split_variables.insert(index);
}

void Data::set_weight_index(size_t index) {
  this->weight_index = index;
  disallowed_split_variables.insert(index);
}

void Data::set_causal_survival_numerator_index(size_t index) {
  this->causal_survival_numerator_index = index;
  disallowed_split_variables.insert(index);
}

void Data::set_causal_survival_denominator_index(size_t index) {
  this->causal_survival_denominator_index = index;
  disallowed_split_variables.insert(index);
}

void Data::set_censor_index(size_t index) {
  this->censor_index = index;
  disallowed_split_variables.insert(index);
}

std::vector<size_t> Data::get_all_values(std::vector<double>& all_values,
                                         std::vector<size_t>& sorted_samples,
                                         const std::vector<size_t>& samples,
                                         size_t var) const {
  all_values.resize(samples.size());
  for (size_t i = 0; i < samples.size(); i++) {
    size_t sample = samples[i];
    all_values[i] = get(sample, var);
  }

  sorted_samples.resize(samples.size());
  std::vector<size_t> index(samples.size());
   // fill with [0, 1,..., samples.size() - 1]
  std::iota(index.begin(), index.end(), 0);
  // sort index based on the split values (argsort)
  // the NaN comparison places all NaNs at the beginning
  // stable sort is needed for consistent element ordering cross platform,
  // otherwise the resulting sums used in the splitting rules may compound rounding error
  // differently and produce different splits.
  std::stable_sort(index.begin(), index.end(), [&](const size_t& lhs, const size_t& rhs) {
    return all_values[lhs] < all_values[rhs] || (std::isnan(all_values[lhs]) && !std::isnan(all_values[rhs]));
  });

  for (size_t i = 0; i < samples.size(); i++) {
    sorted_samples[i] = samples[index[i]];
    all_values[i] = get(sorted_samples[i], var);
  }

  all_values.erase(unique(all_values.begin(), all_values.end(), [&](const double& lhs, const double& rhs) {
    return lhs == rhs || (std::isnan(lhs) && std::isnan(rhs));
  }), all_values.end());

  return index;
}

size_t Data::get_num_cols() const {
  return num_cols;
}

size_t Data::get_num_rows() const {
  return num_rows;
}

size_t Data::get_num_outcomes() const {
  if (outcome_index.has_value()) {
    return outcome_index.value().size();
  } else {
    return 1;
  }
}

size_t Data::get_num_treatments() const {
  if (treatment_index.has_value()) {
    return treatment_index.value().size();
  } else {
    return 1;
  }
}

const std::set<size_t>& Data::get_disallowed_split_variables() const {
  return disallowed_split_variables;
}

} // namespace grf
