/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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

#ifndef GRF_RELABELINGSTRATEGY_H
#define GRF_RELABELINGSTRATEGY_H

#include <vector>

#include "Eigen/Dense"
#include "commons/Data.h"

namespace grf {

/**
 * Produces a relabelled set of outcomes for a set of training samples. These outcomes
 * will then be used in calculating a standard regression (or classification) split.
 *
 * The optional override `get_response_length()` is used for forests splitting on
 * vector-valued outcomes.
 */
class RelabelingStrategy {
public:

  virtual ~RelabelingStrategy() = default;

 /**
   * samples: the subset of samples to relabel.
   * data: the training data matrix.
   * responses_by_sample: the output of the method, an array of relabelled response for each sample ID in `samples`.
   * The dimension of this array is N * K where N is the total number of samples in the data, and K is given
   * by `get_response_length()`.
   *
   * In most cases, like a single-variable regression forest, K is 1, and `responses_by_sample` is a scalar for
   * each sample. In other forests, like multi-output regression forest, K is equal to the number of outcomes,
   * and `responses_by_sample` is a length K vector for each sample (working with a vector-valued splitting rule).
   *
   * Note that for performance reasons (avoiding clearing out the array after each split) this array may
   * contain garbage values for indices outside of the given set of sample IDs.
   *
   * returns: a boolean that will be 'true' if splitting should stop early.
   */
  virtual bool relabel(const std::vector<size_t>& samples,
                       const Data& data,
                       Eigen::ArrayXXd& responses_by_sample) const = 0;

 /**
   * Override to specify the column dimension of `responses_by_sample`.
   * The default value of 1 is used for most forests splitting on scalar values.
   */
  virtual size_t get_response_length() const { return 1; };
};

} // namespace grf

#endif //GRF_RELABELINGSTRATEGY_H
