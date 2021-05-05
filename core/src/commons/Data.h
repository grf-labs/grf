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

#include <set>
#include <vector>

#include "Eigen/Dense"
#include "globals.h"
#include "optional/optional.hpp"

namespace grf {

/**
 * Data wrapper for GRF.
 * Serves as a read-only (immutable) wrapper of a column major (Fortran order)
 * array accessed through its pointer (data_ptr). This class does not own
 * data.
 *
 * The GRF data model is a contiguous array [X, Y, z, ...] of covariates X,
 * outcomes Y, and other optional variables z.
 *
 */
class Data {
public:
  Data(const double* data_ptr, size_t num_rows, size_t num_cols);

  /**
   * Convenience constructors for unit test.
   * The intended use case is with storage (data vector) mananaged
   * elsewhere, i.e.
   * std::vector<double> data_vector {1, 2, 3, etc};
   * Data grf_data_wrapper(data_vector, num_rows, num_cols);
   */
  Data(const std::vector<double>& data, size_t num_rows, size_t num_cols);

  Data(const std::pair<std::vector<double>, std::vector<size_t>>& data);

  void set_outcome_index(size_t index);

  void set_outcome_index(const std::vector<size_t>& index);

  void set_treatment_index(size_t index);

  void set_treatment_index(const std::vector<size_t>& index);

  void set_instrument_index(size_t index);

  void set_weight_index(size_t index);

  void set_causal_survival_numerator_index(size_t index);

  void set_causal_survival_denominator_index(size_t index);

  void set_censor_index(size_t index);

  /**
   * Sorts and gets the unique values in `samples` at variable `var`.
   *
   * @param all_values: the unique values in sorted order (filled in place).
   * @param sorted_samples: the sample IDs in sorted order (filled in place).
   * @param samples: the samples to sort.
   * @param var: the feature variable.
   * @return: (optional) the index (arg sort) of `sorted_samples` (integers from 0,...,samples.size() - 1).
   *
   * If all the values in `samples` is unique, then `all_values` and `sorted_samples`
   * have the same length.
   *
   * If any of the covariates are NaN, they will be placed first in the returned sort order.
   */
  std::vector<size_t> get_all_values(std::vector<double>& all_values,
                                     std::vector<size_t>& sorted_samples,
                                     const std::vector<size_t>& samples, size_t var) const;

  size_t get_num_cols() const;

  size_t get_num_rows() const;

  size_t get_num_outcomes() const;

  size_t get_num_treatments() const;

  const std::set<size_t>& get_disallowed_split_variables() const;

  double get_outcome(size_t row) const;

  Eigen::VectorXd get_outcomes(size_t row) const;

  double get_treatment(size_t row) const;

  Eigen::VectorXd get_treatments(size_t row) const;

  double get_instrument(size_t row) const;

  double get_weight(size_t row) const;

  double get_causal_survival_numerator(size_t row) const;

  double get_causal_survival_denominator(size_t row) const;

  bool is_failure(size_t row) const;

  double get(size_t row, size_t col) const;

private:
  const double* data_ptr;
  size_t num_rows;
  size_t num_cols;

  std::set<size_t> disallowed_split_variables;
  nonstd::optional<std::vector<size_t>> outcome_index;
  nonstd::optional<std::vector<size_t>> treatment_index;
  nonstd::optional<size_t> instrument_index;
  nonstd::optional<size_t> weight_index;
  nonstd::optional<size_t> causal_survival_numerator_index;
  nonstd::optional<size_t> causal_survival_denominator_index;
  nonstd::optional<size_t> censor_index;
};

// inline appropriate getters
inline double Data::get_outcome(size_t row) const {
  return get(row, outcome_index.value()[0]);
}

inline Eigen::VectorXd Data::get_outcomes(size_t row) const {
  Eigen::VectorXd out(outcome_index.value().size());
  for (size_t i = 0; i < outcome_index.value().size(); i++) {
    out(i) = get(row, outcome_index.value()[i]);
  }
  return out;
}

inline double Data::get_treatment(size_t row) const {
  return get(row, treatment_index.value()[0]);
}

inline Eigen::VectorXd Data::get_treatments(size_t row) const {
  Eigen::VectorXd out(treatment_index.value().size());
  for (size_t i = 0; i < treatment_index.value().size(); i++) {
    out(i) = get(row, treatment_index.value()[i]);
  }
  return out;
}

inline double Data::get_instrument(size_t row) const {
  return get(row, instrument_index.value());
}

inline double Data::get_weight(size_t row) const {
  if (weight_index.has_value()) {
    return get(row, weight_index.value());
  } else {
    return 1.0;
  }
}

inline double Data::get_causal_survival_numerator(size_t row) const {
  return get(row, causal_survival_numerator_index.value());
}

inline double Data::get_causal_survival_denominator(size_t row) const {
  return get(row, causal_survival_denominator_index.value());
}

inline bool Data::is_failure(size_t row) const {
  return get(row, censor_index.value()) > 0.0;
}

inline double Data::get(size_t row, size_t col) const {
  return data_ptr[col * num_rows + row];
}

} // namespace grf
#endif /* GRF_DATA_H_ */
