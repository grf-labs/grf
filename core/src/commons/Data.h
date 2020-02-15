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

#ifndef GRF_DATA_H_
#define GRF_DATA_H_

#include <iostream>
#include <set>
#include <vector>

#include "globals.h"
#include "optional/optional.hpp"

namespace grf {

class Data {
public:
  Data();

  virtual ~Data() = default;

  virtual void reserve_memory() = 0;

  virtual double get(size_t row, size_t col) const = 0;

  virtual void set(size_t col, size_t row, double value, bool& error) = 0;

  /**
   * sort() is only needed for the biq Q splitting rule to build an index based on
   * the global sort order. When called, a sweep sets the member `contains_nan()`
   * to true if a NaN is encountered, in which case only the small q splitting should
   * be called because the proceeding sort does not account for NaNs.
   *
   * This means the member functions
   *
   * get_index
   * get_unique_data_value
   * get_num_unique_data_values
   * get_max_num_unique_values
   *
   * can still be used, but should not be relied on for correctness.
   *
   */
  void sort();

  bool load_from_file(const std::string& filename);

  bool load_from_whitespace_file(std::ifstream& input_file, const std::string& first_line);

  bool load_from_other_file(std::ifstream& input_file, const std::string& first_line, char seperator);

  void set_outcome_index(size_t index);

  void set_treatment_index(size_t index);

  void set_instrument_index(size_t index);

  void set_weight_index(size_t index);

  /**
   * Sorts and gets the unique values in `samples` at variable `var`.
   *
   * @param all_values: the unique values in sorted order (filled in place).
   * @param sorted_samples: the sample IDs in sorted order (filled in place).
   * @param samples: the samples to sort.
   * @param var: the feature variable.
   *
   * If all the values in `samples` is unique, then `all_values` and `sorted_samples`
   * have the same length.
   */
  void get_all_values(std::vector<double>& all_values, std::vector<size_t>& sorted_samples, const std::vector<size_t>& samples, size_t var) const;

  size_t get_index(size_t row, size_t col) const;

  double get_unique_data_value(size_t var, size_t index) const;

  size_t get_num_unique_data_values(size_t var) const;

  size_t get_num_cols() const;

  size_t get_num_rows() const;

  size_t get_max_num_unique_values() const;

  double get_outcome(size_t row) const;

  double get_treatment(size_t row) const;

  double get_instrument(size_t row) const;

  double get_weight(size_t row) const;

  const std::set<size_t>& get_disallowed_split_variables() const;

  bool contains_nan() const;

protected:
  size_t num_rows;
  size_t num_cols;
  bool has_nan;

  std::vector<size_t> index_data;
  std::vector<std::vector<double>> unique_data_values;
  size_t max_num_unique_values;

  std::set<size_t> disallowed_split_variables;
  nonstd::optional<size_t> outcome_index;
  nonstd::optional<size_t> treatment_index;
  nonstd::optional<size_t> instrument_index;
  nonstd::optional<size_t> weight_index;

private:
  DISALLOW_COPY_AND_ASSIGN(Data);
};

} // namespace grf
#endif /* GRF_DATA_H_ */
